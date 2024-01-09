/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <fstream>
#include <iostream>
#include <unordered_set>

#include <Eigen/SparseCore>

#include <sym/pose3.h>

#include <lcmtypes/sym/bal_test_data_t.hpp>

#include "gen/snavely_reprojection_factor.h"

#include "alloc.h"
#include "arena.h"
#include "linearizer.h"

void CheckProblem(sym::bal_test_data_t ref, const std::string& filename, std::vector<int> observation_indices) {
  std::ifstream file(filename);

  int num_cameras, num_points, num_observations;
  file >> num_cameras;
  file >> num_points;
  file >> num_observations;

  if (observation_indices.empty()) {
    for (int i = 0; i < num_observations; i++) {
      observation_indices.push_back(i);
    }
  }

  std::vector<int> camera_indices;
  std::vector<int> point_indices;
  std::vector<Eigen::Vector2d> pixels;

  for (int i = 0; i < num_observations; i++) {
    int camera, point;
    file >> camera;
    file >> point;

    double px, py;
    file >> px;
    file >> py;

    camera_indices.push_back(camera);
    point_indices.push_back(point);
    pixels.emplace_back(px, py);
  }

  std::vector<sym::Pose3d> camera_poses;
  std::vector<Eigen::Vector3d> intrinsics;

  for (int i = 0; i < num_cameras; i++) {
    double rx, ry, rz, tx, ty, tz, f, k1, k2;
    file >> rx;
    file >> ry;
    file >> rz;
    file >> tx;
    file >> ty;
    file >> tz;
    file >> f;
    file >> k1;
    file >> k2;

    camera_poses.emplace_back(sym::Rot3d::FromTangent(Eigen::Vector3d(rx, ry, rz)),
                              Eigen::Vector3d(tx, ty, tz));
    intrinsics.emplace_back(f, k1, k2);
  }

  std::vector<Eigen::Vector3d> points;

  for (int i = 0; i < num_points; i++) {
    double x, y, z;
    file >> x;
    file >> y;
    file >> z;

    points.emplace_back(x, y, z);
  }

  // Extract a subset of the full problem and remap the indices.
  std::unordered_map<i32, i32> camera_indices_map;
  for (const auto obs_index : observation_indices) {
    if (camera_indices_map.find(camera_indices.at(obs_index)) == camera_indices_map.end()) {
      camera_indices_map[camera_indices.at(obs_index)] = camera_indices_map.size();
    }
  }

  std::unordered_set<i32> observation_indices_set(observation_indices.begin(), observation_indices.end());

  std::unordered_map<i32, i32> point_indices_map;
  for (const auto obs_index : observation_indices) {
    if (point_indices_map.find(point_indices.at(obs_index)) == point_indices_map.end()) {
      point_indices_map[point_indices.at(obs_index)] = point_indices_map.size();
    }
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
  i32 nblocks = observation_indices.size() * 6;
  i32* Hl_block_rows = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  i32* Hl_block_cols = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);

  // Linearizer callers are responsible for choosing the order of the keys.
  // Here we have all the cameras in order (pose then intrinsics) followed by all the points.
  i32 nkeys = 2 * camera_indices_map.size() + point_indices_map.size();
  for (i32 factor_index = 0; factor_index < observation_indices.size(); ++factor_index) {
    i32 obs_index = observation_indices.at(factor_index);
    i32 camera_index = camera_indices_map.at(camera_indices.at(obs_index));
    i32 pose_key = 2 * camera_index + 0;
    i32 intrinsics_key = 2 * camera_index + 1;
    i32 point_index = point_indices_map.at(point_indices.at(obs_index));
    i32 point_key = 2 * camera_indices_map.size() + point_index;

    // NOTE: The order here must be consistent with the order later on in the update calls.
    Hl_block_rows[6 * factor_index + 0] = pose_key;
    Hl_block_rows[6 * factor_index + 1] = intrinsics_key;
    Hl_block_rows[6 * factor_index + 2] = point_key;
    Hl_block_rows[6 * factor_index + 3] = intrinsics_key;
    Hl_block_rows[6 * factor_index + 4] = point_key;
    Hl_block_rows[6 * factor_index + 5] = point_key;

    Hl_block_cols[6 * factor_index + 0] = pose_key;
    Hl_block_cols[6 * factor_index + 1] = pose_key;
    Hl_block_cols[6 * factor_index + 2] = pose_key;
    Hl_block_cols[6 * factor_index + 3] = intrinsics_key;
    Hl_block_cols[6 * factor_index + 4] = intrinsics_key;
    Hl_block_cols[6 * factor_index + 5] = point_key;
  }

  // Compute key sizes.
  i32* key_sizes = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
  for (i32 i = 0; i < camera_indices_map.size(); ++i) {
    key_sizes[2 * i + 0] = 6;
    key_sizes[2 * i + 1] = 3;
  }
  for (i32 i = 0; i < point_indices_map.size(); ++i) {
    key_sizes[2 * camera_indices_map.size() + i] = 3;
  }

  // Create the linearizer, linearization.
  i32* Hl_block_nz_indices = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
  sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);
  alloc->free(Hl_block_rows, nblocks * sizeof(i32), alloc->ctx);
  alloc->free(Hl_block_cols, nblocks * sizeof(i32), alloc->ctx);

  i32* key_perm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
  i32 use_metis = 1;
  if (use_metis) {
      sym_get_metis_lower_tri_perm(Hl_block, key_sizes, key_perm, alloc);
  } else {
      for (i32 i = 0; i < nkeys; ++i) {
          key_perm[i] = i;
      }
  }

  i32* key_size_scan = (i32*) alloc->malloc((nkeys + 1) * sizeof(i32), alloc->ctx);
  key_size_scan[0] = 0;
  for (i32 i = 0; i < nkeys; ++i) {
      key_size_scan[i + 1] = key_size_scan[i] + key_sizes[i];
  }

  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,i32> P(key_size_scan[nkeys]);
  i32 pindex = 0;
  for (i32 i = 0; i < nkeys; ++i) {
      i32 j = key_perm[i];
      i32 key_size = key_sizes[j];
      i32 offset = key_size_scan[j];
      for (i32 k = 0; k < key_size; ++k) {
          P.indices()(offset + k) = pindex++;
      }
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

  // Run a single linearization round.
  sym_linearization_clear(lin);

  // (cam_T_world: 6, intrinsics: 3, point: 3)
  double fac_res[2];
  double fac_hessian_dense[12 * 12];
  double fac_rhs[12];

  f64 error = 0.0;
  for (i32 factor_index = 0; factor_index < observation_indices.size(); factor_index++) {
    i32 obs_index = observation_indices.at(factor_index);
    i32 camera_index = camera_indices_map.at(camera_indices.at(obs_index));
    i32 pose_key = 2 * camera_index + 0;
    i32 intrinsics_key = 2 * camera_index + 1;
    i32 point_index = point_indices_map.at(point_indices.at(obs_index));
    i32 point_key = 2 * camera_indices_map.size() + point_index;

    sym::Pose3d camera_pose = camera_poses.at(camera_indices.at(obs_index));
    Eigen::Vector3d intrinsics_vec = intrinsics.at(camera_indices.at(obs_index));
    Eigen::Vector3d point = points.at(point_indices.at(obs_index));
    Eigen::Vector2d pixels_vec = pixels.at(obs_index);

    double camera_pose_storage[7];
    camera_pose.ToStorage(camera_pose_storage);

    snavely_reprojection_factor(camera_pose_storage, 
                              intrinsics_vec.data(),
                              point.data(),
                              pixels_vec.data(),
                              sym::kDefaultEpsilond,
                              // is there a better way to pass in nullptr?
                              fac_res, NULL, fac_hessian_dense, fac_rhs);
                            
    error += fac_res[0] * fac_res[0];
    error += fac_res[1] * fac_res[1];

    // NOTE: The block index order must match the order in the block triplets above.
    sym_linearizer_add_hessian_tri_block(
      lzr, lin, 
      6 * factor_index + 0, pose_key, 
      fac_hessian_dense, 12, 0
    );
    sym_linearizer_add_hessian_rect_block(
      lzr, lin, 
      6 * factor_index + 1, intrinsics_key, pose_key,
      fac_hessian_dense, 12, 6, 0
    );
    sym_linearizer_add_hessian_rect_block(
      lzr, lin, 
      6 * factor_index + 2, point_key, pose_key,
      fac_hessian_dense, 12, 9, 0
    );

    sym_linearizer_add_hessian_tri_block(
      lzr, lin, 
      6 * factor_index + 3, intrinsics_key, 
      fac_hessian_dense, 12, 6
    );
    sym_linearizer_add_hessian_rect_block(
      lzr, lin, 
      6 * factor_index + 4, point_key, intrinsics_key,
      fac_hessian_dense, 12, 9, 6
    );

    sym_linearizer_add_hessian_tri_block(
      lzr, lin, 
      6 * factor_index + 5, point_key, 
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
  error *= 0.5;

  Eigen::Map<Eigen::SparseMatrix<double>> Hl(
    lin.Hl.nrows,
    lin.Hl.ncols,
    lin.Hl.nnz,
    lin.Hl.col_starts,
    lin.Hl.row_indices,
    lin.Hl.data
  );
  Eigen::Map<Eigen::VectorXd> rhs(
    lin.rhs.data,
    lin.rhs.n
  );

  // Load the reference Hessian.
  // NOTE: a nnz entry needs to be added to the column pointers array
  std::vector<int> column_pointers(ref.hessian_structure.shape[1] + 1);
  std::memcpy(column_pointers.data(), ref.hessian_structure.column_pointers.data(), column_pointers.size() * sizeof(int));
  column_pointers.back() = ref.hessian_structure.row_indices.size();
  Eigen::Map<const Eigen::SparseMatrix<double>> Hl_ref(
    ref.hessian_structure.shape[0],
    ref.hessian_structure.shape[1],
    ref.hessian_structure.row_indices.size(),
    column_pointers.data(),
    ref.hessian_structure.row_indices.data(),
    ref.hessian_data.data()
  );

  Eigen::SparseMatrix<f64> Hl_ref_perm;
  Hl_ref_perm.template selfadjointView<Eigen::Lower>() =
        Hl_ref.template selfadjointView<Eigen::Lower>().twistedBy(P);
  // TODO: is the easiest way to sort the entries in each col?
  Hl_ref_perm = Hl_ref_perm.transpose();
  Hl_ref_perm = Hl_ref_perm.transpose();

  // Compare!
  assert(Hl.rows() == Hl_ref.rows());
  assert(Hl.cols() == Hl_ref.cols());
  assert(Hl.nonZeros() == Hl_ref.nonZeros());

  assert(Hl.nonZeros() == Hl_ref_perm.nonZeros());

  if (!Hl.isApprox(Hl_ref_perm)) {
    std::cout <<  "hessian: " << std::endl << Hl << std::endl << std::endl;
    std::cout << "ref_hessian: " << std::endl << Eigen::Map<const Eigen::SparseMatrix<double>>(Hl_ref_perm.rows(), 
      Hl_ref_perm.cols(), Hl_ref_perm.nonZeros(), Hl_ref_perm.outerIndexPtr(), Hl_ref_perm.innerIndexPtr(), Hl_ref_perm.valuePtr()) << std::endl << std::endl;
    std::cout << "hessian - ref_hessian: " << (Hl - Hl_ref_perm) << std::endl;
    assert(false);
  }

  Eigen::VectorXd rhs_ref_perm = P * ref.rhs;

  assert(rhs.size() == ref.rhs.size());
  if (!rhs.isApprox(rhs_ref_perm)) {
    std::cout << "rhs: " << rhs << std::endl;
    std::cout << "ref_rhs: " << rhs_ref_perm << std::endl;
    std::cout << "rhs - ref_rhs: " << (rhs - rhs_ref_perm) << std::endl;
    assert(false);
  }

  sym_linearizer_free(lzr, alloc);
  sym_linearization_free(lin, alloc);

  alloc->free(Hl_block_nz_indices, nblocks * sizeof(i32), alloc->ctx);

  printf("nalloc = %d\n", arena.nalloc);
  printf("max_nalloc = %d\n", arena.max_nalloc);
  assert(arena.nalloc == 0);
  
  free(buf);
}

sym::bal_test_data_t ReadResults(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  std::vector<uint8_t> buffer(std::istreambuf_iterator<char>(file), {});
  sym::bal_test_data_t data;
  data.decode(buffer.data(), 0, buffer.size());
  return data;
}

int main(int argc, char** argv) {
  const auto* filename = "test/bal/data/problem-21-11315-pre.txt";
  CheckProblem(ReadResults("test/bal/data_ref/0.bin"), filename, {0});
  CheckProblem(ReadResults("test/bal/data_ref/1.bin"), filename, {0, 0});
  CheckProblem(ReadResults("test/bal/data_ref/2.bin"), filename, {0, 1});
  CheckProblem(ReadResults("test/bal/data_ref/3.bin"), filename, {13, 14});
  CheckProblem(ReadResults("test/bal/data_ref/4.bin"), filename, {});
}
