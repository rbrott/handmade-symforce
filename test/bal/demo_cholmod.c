/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <math.h>
#include <float.h>
#include <stdbool.h>

#include <cholmod.h>

#include "gen/snavely_reprojection_factor.h"
#include "gen/pose3_retract.h"
#include "gen/rot3_tangent.h"

#include "sym_assert.h"
#include "alloc.h"
#include "arena.h"
#include "linearizer.h"

typedef struct {
  i32* camera_indices;
  i32* point_indices;
  f64* pixels;
  i32 num_cameras;
  i32 num_points;
  i32 num_observations;
} bal_problem;

void print_error(int status, const char *file, int line,
      const char *message) {
  printf("status: %d, file: %s, line: %d, message: %s\n", status, file, line, message);
}

f64 bal_linearize(bal_problem p, sym_linearizer lzr, sym_linearization lin, f64* values, f64 epsilon) {
  // Run a single linearization round.
  sym_linearization_clear(lin);

  // (cam_T_world: 6, intrinsics: 3, point: 3)
  f64 fac_res[2];
  f64 fac_hessian_dense[12 * 12];
  f64 fac_rhs[12];

  f64 error = 0.0;
  for (i32 obs_index = 0; obs_index < p.num_observations; obs_index++) {
    i32 camera_index = p.camera_indices[obs_index];
    i32 pose_key = 2 * camera_index + 0;
    i32 intrinsics_key = 2 * camera_index + 1;
    i32 point_index = p.point_indices[obs_index];
    i32 point_key = 2 * p.num_cameras + point_index;

    snavely_reprojection_factor(values + 10 * camera_index,
                                values + 10 * camera_index + 7,
                                values + 10 * p.num_cameras + 3 * point_index,
                                p.pixels + 2 * obs_index,
                                epsilon,
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
}

int main(int argc, char** argv) {
  SYM_ASSERT(argc == 2);

  f64 epsilon = 10.0 * DBL_EPSILON;

  size n = 1L << 30; // 1 GiB seems fine
  printf("Allocating %ld bytes\n", n);
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

  FILE* file = fopen(argv[1], "r");

  bal_problem p = {};
  fscanf(file, "%d", &p.num_cameras);
  fscanf(file, "%d", &p.num_points);
  fscanf(file, "%d", &p.num_observations);

  p.camera_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);
  p.point_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);
  // x, y interleaved
  p.pixels = (f64*) alloc->malloc(2 * p.num_observations * sizeof(f64), alloc->ctx);

  for (i32 i = 0; i < p.num_observations; i++) {
    i32 camera, point;
    fscanf(file, "%d", &camera);
    fscanf(file, "%d", &point);

    f64 px, py;
    fscanf(file, "%lf", &px);
    fscanf(file, "%lf", &py);

    p.camera_indices[i] = camera;
    p.point_indices[i] = point;
    p.pixels[2 * i + 0] = px;
    p.pixels[2 * i + 1] = py;
  }

  // (cam_T_world: 7, intrinsics: 3, point: 3)
  i32 values_dim = (7 + 3) * p.num_cameras + 3 * p.num_points;
  f64* values = (f64*) alloc->malloc(values_dim * sizeof(f64), alloc->ctx);
  f64* temp_values = (f64*) alloc->malloc(values_dim * sizeof(f64), alloc->ctx);
  for (i32 i = 0; i < p.num_cameras; i++) {
    f64 rx, ry, rz, tx, ty, tz, f, k1, k2;
    fscanf(file, "%lf", &rx);
    fscanf(file, "%lf", &ry);
    fscanf(file, "%lf", &rz);
    fscanf(file, "%lf", &tx);
    fscanf(file, "%lf", &ty);
    fscanf(file, "%lf", &tz);
    fscanf(file, "%lf", &f);
    fscanf(file, "%lf", &k1);
    fscanf(file, "%lf", &k2);

    f64 rot_tangent[3] = {rx, ry, rz};
    rot3_from_tangent(rot_tangent, values + 10 * i, epsilon);

    values[10 * i + 4] = tx;
    values[10 * i + 5] = ty;
    values[10 * i + 6] = tz;

    values[10 * i + 7] = f;
    values[10 * i + 8] = k1;
    values[10 * i + 9] = k2;
  }

  for (i32 i = 0; i < p.num_points; i++) {
    f64 x, y, z;
    fscanf(file, "%lf", &x);
    fscanf(file, "%lf", &y);
    fscanf(file, "%lf", &z);

    values[10 * p.num_cameras + 3 * i + 0] = x;
    values[10 * p.num_cameras + 3 * i + 1] = y;
    values[10 * p.num_cameras + 3 * i + 2] = z;
  }

  fclose(file);

  // Compute Hessian_lower block triplets.
  i32 nblocks = p.num_observations * 6;
  i32* Hl_block_rows = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  i32* Hl_block_cols = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);

  // Linearizer callers are responsible for choosing the order of the keys.
  // Here we have all the cameras in order (pose then intrinsics) followed by all the points.
  i32 nkeys = 2 * p.num_cameras + p.num_points;
  for (i32 obs_index = 0; obs_index < p.num_observations; ++obs_index) {
    i32 camera_index = p.camera_indices[obs_index];
    i32 pose_key = 2 * camera_index + 0;
    i32 intrinsics_key = 2 * camera_index + 1;
    i32 point_index = p.point_indices[obs_index];
    i32 point_key = 2 * p.num_cameras + point_index;

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
  for (i32 i = 0; i < p.num_cameras; ++i) {
    key_sizes[2 * i + 0] = 6;
    key_sizes[2 * i + 1] = 3;
  }
  for (i32 i = 0; i < p.num_points; ++i) {
    key_sizes[2 * p.num_cameras + i] = 3;
  }

  // Create the linearizer, linearization.
  i32* Hl_block_nz_indices = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);
  alloc->free(Hl_block_rows, nblocks * sizeof(i32), alloc->ctx);
  alloc->free(Hl_block_cols, nblocks * sizeof(i32), alloc->ctx);

  i32* key_perm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
  sym_get_metis_tri_perm(Hl_block, key_sizes, key_perm, alloc);

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

  cholmod_common chol_common;
  cholmod_start(&chol_common);
  chol_common.nmethods = 1;
  chol_common.method[0].ordering = CHOLMOD_NATURAL;

  chol_common.error_handler = print_error;


  cholmod_sparse H = {
    .nrow = lin.Hl.nrows,
    .ncol = lin.Hl.ncols,
    .nzmax = lin.Hl.nnz,
    .p = lin.Hl.col_starts,
    .i = lin.Hl.row_indices,
    .nz = NULL,
    .x = NULL,
    .z = NULL,
    .stype = -1, // symmetric, lower part stored
    .itype = CHOLMOD_INT,
    .xtype = CHOLMOD_PATTERN,
    .dtype = CHOLMOD_DOUBLE,
    .sorted = 1,
    .packed = 1,
  };

  // H, H_fac symbolic for now -- does this allocate memory for the data members?
  cholmod_factor* H_fac = cholmod_analyze(&H, &chol_common);

// typedef struct cholmod_dense_struct
// {
//     size_t nrow ;       // the matrix is nrow-by-ncol
//     size_t ncol ;
//     size_t nzmax ;      // maximum number of entries in the matrix
//     size_t d ;          // leading dimension (d >= nrow must hold)
//     void *x ;           // size nzmax or 2*nzmax, if present
//     void *z ;           // size nzmax, if present
//     int xtype ;         // pattern, real, complex, or zomplex
//     int dtype ;         // x and z double or single

// } cholmod_dense ;

  cholmod_dense b = {
      .nrow = lin.Hl.nrows,
      .ncol = 1,
      .nzmax = lin.Hl.nrows,
      .d = lin.Hl.nrows,
      .x = alloc->malloc(lin.Hl.nrows * sizeof(f64), alloc->ctx),
      .z = NULL,
      .xtype = CHOLMOD_REAL,
      .dtype = CHOLMOD_DOUBLE
  };

  // LM params
  f64 initial_lambda = 1.0;
  f64 lambda_up_factor = 4.0;
  f64 lambda_down_factor = 1 / 4.0;
  f64 lambda_lower_bound = 0.0;
  f64 lambda_upper_bound = 1000000.0;
  f64 early_exit_min_reduction = 1e-6;

  f64 last_error = bal_linearize(p, lzr, lin, values, epsilon);
  f64 lambda = initial_lambda;
  i32 iteration = 0;
  while (1) {
    // damp the last Hessian (just unit diagonal damping)
    for (i32 i = 0; i < lin.Hl.nrows; ++i) {
      i32 j = lin.Hl.col_starts[i];
      SYM_ASSERT(lin.Hl.row_indices[j] == i);
      lin.Hl.data[j] += lambda;
    }

    H.xtype = CHOLMOD_REAL;
    H.x = lin.Hl.data;

    // TODO: use the fancier routines to avoid allocating (like solve2)

    int success = cholmod_factorize(&H, H_fac, &chol_common);
    SYM_ASSERT(success); // TODO: idk what the value of this is supposed to be

    for (i32 i = 0; i < lin.Hl.nrows; ++i) {
      ((f64*) b.x)[i] = -lin.rhs.data[i];
    }

    cholmod_dense* x = cholmod_solve(CHOLMOD_A, H_fac, &b, &chol_common);

    // copy values into temp
    for (i32 i = 0; i < values_dim; ++i) {
      temp_values[i] = values[i];
    }

    // apply the update to the temp values
    f64* xdata = (f64*) x->x;
    for (i32 i = 0; i < p.num_cameras; ++i) {
      i32 pose_key = 2 * i + 0;
      i32 pose_values_offset = 10 * i;
      i32 pose_rhs_offset = lzr.key_size_scan[lzr.key_iperm[pose_key]];
      sym_pose3_retract_in_place(temp_values + pose_values_offset, xdata + pose_rhs_offset, epsilon);

      i32 intrinsics_key = 2 * i + 1;
      i32 intrinsics_values_offset = 10 * i + 7;
      i32 intrinsics_rhs_offset = lzr.key_size_scan[lzr.key_iperm[intrinsics_key]];
      for (i32 j = 0; j < 3; ++j) {
        temp_values[intrinsics_values_offset + j] += xdata[intrinsics_rhs_offset + j];
      }
    }
    for (i32 i = 0; i < p.num_points; ++i) {
      i32 point_key = 2 * p.num_cameras + i;
      i32 point_values_offset = 10 * p.num_cameras + 3 * i;
      i32 point_rhs_offset = lzr.key_size_scan[lzr.key_iperm[point_key]];
      for (i32 j = 0; j < 3; ++j) {
        temp_values[point_values_offset + j] += xdata[point_rhs_offset + j];
      }
    }

    f64 error = bal_linearize(p, lzr, lin, temp_values, epsilon);
    f64 relative_reduction = (last_error - error) / (last_error + epsilon);

    printf("BAL optimizer [iter %4d] lambda: %e, error prev/new: %e/%e, rel reduction: %e\n", 
      iteration, lambda, last_error, error, relative_reduction);

    if (relative_reduction > -early_exit_min_reduction / 10 &&
        relative_reduction < early_exit_min_reduction) {
      // TODO: do we care about setting the values?
      printf("Success!\n");
      break;
    }

    bool accept_update = relative_reduction > 0;

    if (!accept_update && lambda >= lambda_upper_bound) {
      printf("Failed: lambda out of bounds!\n");
      break;
    }

    if (!accept_update) {
      lambda *= lambda_up_factor;
    } else {
      lambda *= lambda_down_factor;
      // TODO: is this right?
      last_error = error;

      // swap values and temp_values
      {
        f64* values_tmp = values;
        values = temp_values;
        temp_values = values_tmp;
      }
    }

    lambda = fmax(fmin(lambda, lambda_upper_bound), lambda_lower_bound);

    ++iteration;
  }

  alloc->free(b.x, lin.Hl.nrows * sizeof(f64), alloc->ctx);

  sym_linearizer_free(lzr, alloc);
  sym_linearization_free(lin, alloc);

  alloc->free(Hl_block_nz_indices, nblocks * sizeof(i32), alloc->ctx);

  alloc->free(p.camera_indices, p.num_observations * sizeof(i32), alloc->ctx);
  alloc->free(p.point_indices, p.num_observations * sizeof(i32), alloc->ctx);
  alloc->free(p.pixels, 2 * p.num_observations * sizeof(f64), alloc->ctx);

  alloc->free(values, values_dim * sizeof(f64), alloc->ctx);
  alloc->free(temp_values, values_dim * sizeof(f64), alloc->ctx);

  printf("nalloc = %d\n", arena.nalloc);
  printf("max_nalloc = %d\n", arena.max_nalloc);
  SYM_ASSERT(arena.nalloc == 0);
  
  free(buf);
}
