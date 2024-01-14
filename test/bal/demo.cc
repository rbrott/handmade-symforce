/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <fstream>
#include <iostream>
#include <unordered_set>

#include <Eigen/Core>

#include <sym/pose3.h>

#include "gen/snavely_reprojection_factor.h"

#include "alloc.h"
#include "arena.h"
#include "linearizer.h"
#include "solver.h"

void RunProblem(const std::string& filename) {
  std::ifstream file(filename);

  int num_cameras, num_points, num_observations;
  file >> num_cameras;
  file >> num_points;
  file >> num_observations;

  i32* camera_indices = (i32*) malloc(num_observations * sizeof(i32));
  i32* point_indices = (i32*) malloc(num_observations * sizeof(i32));
  // x, y interleaved
  f64* pixels = (f64*) malloc(2 * num_observations * sizeof(f64));

  for (i32 i = 0; i < num_observations; i++) {
    i32 camera, point;
    file >> camera;
    file >> point;

    f64 px, py;
    file >> px;
    file >> py;

    camera_indices[i] = camera;
    point_indices[i] = point;
    pixels[2 * i + 0] = px;
    pixels[2 * i + 1] = py;
  }

  // (cam_T_world: 7, intrinsics: 3, point: 3)
  f64* values = (f64*) malloc(((7 + 3) * num_cameras + 3 * num_points) * sizeof(f64));
  for (i32 i = 0; i < num_cameras; i++) {
    f64 rx, ry, rz, tx, ty, tz, f, k1, k2;
    file >> rx;
    file >> ry;
    file >> rz;
    file >> tx;
    file >> ty;
    file >> tz;
    file >> f;
    file >> k1;
    file >> k2;

    sym::Pose3d cam_T_world(sym::Rot3d::FromTangent(Eigen::Vector3d(rx, ry, rz)),
                              Eigen::Vector3d(tx, ty, tz));
    cam_T_world.ToStorage(values + 10 * i);

    values[10 * i + 7] = f;
    values[10 * i + 8] = k1;
    values[10 * i + 9] = k2;
  }

  for (i32 i = 0; i < num_points; i++) {
    f64 x, y, z;
    file >> x;
    file >> y;
    file >> z;

    values[10 * num_cameras + 3 * i + 0] = x;
    values[10 * num_cameras + 3 * i + 1] = y;
    values[10 * num_cameras + 3 * i + 2] = z;
  }

  size n = 1 << 30; // 1 GiB seems fine
  u8* buf = (u8*) malloc(n);

  sym_arena arena = {
      .beg = buf,
      .end = buf + n,
  };

  sym_allocator alloc_struct = {
      .malloc = sym_arena_malloc,
      .free = sym_arena_free,
      .ctx = &arena
  };
  sym_allocator* alloc = &alloc_struct;

  // Compute Hessian_lower block triplets.
  i32 nblocks = num_observations * 6;
  i32* Hl_block_rows = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  i32* Hl_block_cols = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);

  // Linearizer callers are responsible for choosing the order of the keys.
  // Here we have all the cameras in order (pose then intrinsics) followed by all the points.
  i32 nkeys = 2 * num_cameras + num_points;
  for (i32 obs_index = 0; obs_index < num_observations; ++obs_index) {
    i32 camera_index = camera_indices[obs_index];
    i32 pose_key = 2 * camera_index + 0;
    i32 intrinsics_key = 2 * camera_index + 1;
    i32 point_index = point_indices[obs_index];
    i32 point_key = 2 * num_cameras + point_index;

    // NOTE: The order here must be consistent with the order later on in the update calls.
    Hl_block_rows[6 * obs_index + 0] = pose_key;
    Hl_block_rows[6 * obs_index + 1] = intrinsics_key;
    Hl_block_rows[6 * obs_index + 2] = point_key;
    Hl_block_rows[6 * obs_index + 3] = intrinsics_key;
    Hl_block_rows[6 * obs_index + 4] = point_key;
    Hl_block_rows[6 * obs_index + 5] = point_key;

    Hl_block_cols[6 * obs_index + 0] = pose_key;
    Hl_block_cols[6 * obs_index + 1] = pose_key;
    Hl_block_cols[6 * obs_index + 2] = pose_key;
    Hl_block_cols[6 * obs_index + 3] = intrinsics_key;
    Hl_block_cols[6 * obs_index + 4] = intrinsics_key;
    Hl_block_cols[6 * obs_index + 5] = point_key;
  }

  // Compute key sizes.
  i32* key_sizes = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
  for (i32 i = 0; i < num_cameras; ++i) {
    key_sizes[2 * i + 0] = 6;
    key_sizes[2 * i + 1] = 3;
  }
  for (i32 i = 0; i < num_points; ++i) {
    key_sizes[2 * num_cameras + i] = 3;
  }

  // Create the linearizer, linearization.
  i32* Hl_block_nz_indices = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);
  alloc->free(Hl_block_rows, nblocks * sizeof(i32), alloc->ctx);
  alloc->free(Hl_block_cols, nblocks * sizeof(i32), alloc->ctx);

  i32* key_perm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
  sym_get_metis_tri_perm(Hl_block, key_sizes, key_perm, alloc);

  i32* key_size_scan = (i32*) alloc->malloc((nkeys + 1) * sizeof(i32), alloc->ctx);
  key_size_scan[0] = 0;
  for (i32 i = 0; i < nkeys; ++i) {
      key_size_scan[i + 1] = key_size_scan[i] + key_sizes[i];
  }

  alloc->free(key_size_scan, (nkeys + 1) * sizeof(i32), alloc->ctx);

  sym_linearization lin;
  sym_linearizer lzr = sym_linearizer_new(
    Hl_block, Hl_block_nz_indices, nblocks,
    key_sizes, nkeys, 
    key_perm,
    &lin,
    alloc
  );

  sym_csc_mat_free(Hl_block, alloc);

  alloc->free(key_perm, nkeys * sizeof(i32), alloc->ctx);
  alloc->free(key_sizes, nkeys * sizeof(i32), alloc->ctx);

  i32* Hlt_perm = (i32*) alloc->malloc(lin.Hl.nnz * sizeof(i32), alloc->ctx);
  sym_csc_mat Hlt = sym_transpose_csc(lin.Hl, Hlt_perm, alloc);

  sym_chol_factorization fac{};
  sym_chol_solver solver = sym_new_chol_solver(Hlt, &fac, alloc);

  sym_vec x = sym_vec_new(lin.Hl.nrows, alloc);

  const auto linearize = [&]() {
    // Run a single linearization round.
    sym_linearization_clear(lin);

    // (cam_T_world: 6, intrinsics: 3, point: 3)
    double fac_res[2];
    double fac_hessian_dense[12 * 12];
    double fac_rhs[12];

    f64 error = 0.0;
    for (i32 obs_index = 0; obs_index < num_observations; obs_index++) {
      i32 camera_index = camera_indices[obs_index];
      i32 pose_key = 2 * camera_index + 0;
      i32 intrinsics_key = 2 * camera_index + 1;
      i32 point_index = point_indices[obs_index];
      i32 point_key = 2 * num_cameras + point_index;

      snavely_reprojection_factor(values + 10 * camera_index,
                                  values + 10 * camera_index + 7,
                                  values + 10 * num_cameras + 3 * point_index,
                                  pixels + 2 * obs_index,
                                  sym::kDefaultEpsilond,
                                  fac_res, NULL, fac_hessian_dense, fac_rhs);
                              
      error += fac_res[0] * fac_res[0];
      error += fac_res[1] * fac_res[1];

      // NOTE: The block index order must match the order in the block triplets above.
      sym_linearizer_add_hessian_tri_block(
        lzr, lin, 
        6 * obs_index + 0, pose_key, 
        fac_hessian_dense, 12, 0
      );
      sym_linearizer_add_hessian_rect_block(
        lzr, lin, 
        6 * obs_index + 1, intrinsics_key, pose_key,
        fac_hessian_dense, 12, 6, 0
      );
      sym_linearizer_add_hessian_rect_block(
        lzr, lin, 
        6 * obs_index + 2, point_key, pose_key,
        fac_hessian_dense, 12, 9, 0
      );

      sym_linearizer_add_hessian_tri_block(
        lzr, lin, 
        6 * obs_index + 3, intrinsics_key, 
        fac_hessian_dense, 12, 6
      );
      sym_linearizer_add_hessian_rect_block(
        lzr, lin, 
        6 * obs_index + 4, point_key, intrinsics_key,
        fac_hessian_dense, 12, 9, 6
      );

      sym_linearizer_add_hessian_tri_block(
        lzr, lin, 
        6 * obs_index + 5, point_key, 
        fac_hessian_dense, 12, 9
      );

      sym_linearizer_add_rhs_block(
        lzr, lin, pose_key, 
        fac_rhs, 0
      );
      sym_linearizer_add_rhs_block(
        lzr, lin, intrinsics_key, 
        fac_rhs, 6
      );
      sym_linearizer_add_rhs_block(
        lzr, lin, point_key, 
        fac_rhs, 9
      );
    }
    return 0.5 * error;
  };

  // default LM params
  const double initial_lambda = 1.0;
  const double lambda_up_factor = 4.0;
  const double lambda_down_factor = 1 / 4.0;
  const double lambda_lower_bound = 0.0;
  const double lambda_upper_bound = 1000000.0;
  const bool use_diagonal_damping = false;
  const bool use_unit_damping = true;
  const bool keep_max_diagonal_damping = false;
  const double diagonal_damping_min = 1e-6;
  const int iterations = 50;
  const double early_exit_min_reduction = 1e-6;
  const bool enable_bold_updates = false;

  f64 last_error = linearize();
  printf("error = %e\n", last_error);
  f64 lambda = initial_lambda;
  i32 iteration = 0;
  while (true) {
    // damp the last Hessian
    assert(use_unit_damping);
    assert(!use_diagonal_damping);
    for (i32 i = 0; i < nkeys; ++i) {
      i32 j = lin.Hl.col_starts[i];
      assert(lin.Hl.row_indices[j] == i);
      lin.Hl.data[j] += lambda;
    }

    // transpose it
    for (i32 i = 0; i < lin.Hl.nnz; ++i) {
      Hlt.data[i] = lin.Hl.data[Hlt_perm[i]];
    }

    sym_chol_solver_factor(solver, Hlt, fac);

    for (i32 i = 0; i < nkeys; ++i) {
      x.data[i] = -lin.rhs.data[i];
    }

    sym_chol_solver_solve_in_place(fac, x, alloc);

    // apply the update to the values

    // printf("error = %e\n", error);

    break;

    ++iteration;
  }

  alloc->free(Hlt_perm, lin.Hl.nnz * sizeof(i32), alloc->ctx);

  sym_vec_free(x, alloc);

  sym_chol_solver_free(solver, alloc);
  sym_chol_factorization_free(fac, alloc);
  sym_csc_mat_free(Hlt, alloc);

  sym_linearizer_free(lzr, alloc);
  sym_linearization_free(lin, alloc);

  alloc->free(Hl_block_nz_indices, nblocks * sizeof(i32), alloc->ctx);

  printf("nalloc = %d\n", arena.nalloc);
  printf("max_nalloc = %d\n", arena.max_nalloc);
  assert(arena.nalloc == 0);
  
  free(buf);
}

int main(int argc, char** argv) {
  assert(argc == 2);

  const char* filename = argv[1];
  RunProblem(filename);
}
