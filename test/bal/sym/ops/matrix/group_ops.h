// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     geo_package/ops/matrix/group_ops.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include "../group_ops.h"

namespace sym {

/**
 * C++ GroupOps implementation for matrices.
 */
template <typename ScalarType, int Rows, int Cols>
struct GroupOps<Eigen::Matrix<ScalarType, Rows, Cols>> {
  using Scalar = ScalarType;
  using T = Eigen::Matrix<Scalar, Rows, Cols>;
  static_assert(std::is_floating_point<ScalarType>::value, "");

  static T Identity() {
    return T::Zero();
  }

  static T Inverse(const T& a) {
    return -a;
  }

  static T Compose(const T& a, const T& b) {
    return b + a;
  }

  static T Between(const T& a, const T& b) {
    return b - a;
  }
};

}  // namespace sym

// Explicit instantiation
extern template struct sym::GroupOps<Eigen::Matrix<double, 1, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 2, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 3, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 4, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 5, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 6, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 7, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 8, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 9, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 2, 2>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 3, 3>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 4, 4>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 5, 5>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 6, 6>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 7, 7>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 8, 8>>;
extern template struct sym::GroupOps<Eigen::Matrix<double, 9, 9>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 1, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 2, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 3, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 4, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 5, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 6, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 7, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 8, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 9, 1>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 2, 2>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 3, 3>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 4, 4>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 5, 5>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 6, 6>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 7, 7>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 8, 8>>;
extern template struct sym::GroupOps<Eigen::Matrix<float, 9, 9>>;