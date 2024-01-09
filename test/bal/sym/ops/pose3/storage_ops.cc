// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     ops/CLASS/storage_ops.cc.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#include "./storage_ops.h"

#include <algorithm>
#include <cassert>

#include <Eigen/Dense>

#include <sym/pose3.h>

namespace sym {

template <typename ScalarType>
void StorageOps<Pose3<ScalarType>>::ToStorage(const Pose3<ScalarType>& a, ScalarType* out) {
  assert(out != nullptr);
  std::copy_n(a.Data().data(), a.StorageDim(), out);
}

template <typename ScalarType>
Pose3<ScalarType> StorageOps<Pose3<ScalarType>>::FromStorage(const ScalarType* data) {
  assert(data != nullptr);
  return Pose3<ScalarType>(Eigen::Map<const typename Pose3<ScalarType>::DataVec>(data));
}

}  // namespace sym

// Explicit instantiation
template struct sym::StorageOps<sym::Pose3<double>>;
template struct sym::StorageOps<sym::Pose3<float>>;