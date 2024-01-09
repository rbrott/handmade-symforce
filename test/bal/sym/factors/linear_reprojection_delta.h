// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/linear_camera_cal.h>
#include <sym/pose3.h>

namespace sym {

/**
 * Reprojects the landmark into the target camera and returns the delta from the correspondence to
 * the reprojection.
 *
 * The landmark is specified as a pixel in the source camera and an inverse range; this means the
 * landmark is fixed in the source camera and always has residual 0 there (this 0 residual is not
 * returned, only the residual in the target camera is returned).
 *
 * Args:
 *     source_pose: The pose of the source camera
 *     source_calibration: The source camera calibration
 *     target_pose: The pose of the target camera
 *     target_calibration: The target camera calibration
 *     source_inverse_range: The inverse range of the landmark in the source camera
 *     source_pixel: The location of the landmark in the source camera
 *     target_pixel: The location of the correspondence in the target camera
 *     epsilon: Small positive value
 *     camera_model_class: The subclass of CameraCal to use as the camera model
 *
 * Outputs:
 *     res: 2dof pixel reprojection error
 *     valid: is valid projection or not
 */
template <typename Scalar>
void LinearReprojectionDelta(const sym::Pose3<Scalar>& source_pose,
                             const sym::LinearCameraCal<Scalar>& source_calibration,
                             const sym::Pose3<Scalar>& target_pose,
                             const sym::LinearCameraCal<Scalar>& target_calibration,
                             const Scalar source_inverse_range,
                             const Eigen::Matrix<Scalar, 2, 1>& source_pixel,
                             const Eigen::Matrix<Scalar, 2, 1>& target_pixel, const Scalar epsilon,
                             Eigen::Matrix<Scalar, 2, 1>* const reprojection_delta = nullptr,
                             Scalar* const is_valid = nullptr) {
  // Total ops: 123

  // Input arrays
  const Eigen::Matrix<Scalar, 7, 1>& _source_pose = source_pose.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _source_calibration = source_calibration.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _target_pose = target_pose.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _target_calibration = target_calibration.Data();

  // Intermediate terms (32)
  const Scalar _tmp0 = 2 * _target_pose[2];
  const Scalar _tmp1 = _target_pose[3] * _tmp0;
  const Scalar _tmp2 = 2 * _target_pose[0] * _target_pose[1];
  const Scalar _tmp3 = 2 * _source_pose[1];
  const Scalar _tmp4 = _source_pose[0] * _tmp3;
  const Scalar _tmp5 = 2 * _source_pose[3];
  const Scalar _tmp6 = _source_pose[2] * _tmp5;
  const Scalar _tmp7 = -_source_calibration[2] + source_pixel(0, 0);
  const Scalar _tmp8 = -_source_calibration[3] + source_pixel(1, 0);
  const Scalar _tmp9 =
      std::pow(Scalar(epsilon + 1 +
                      std::pow(_tmp8, Scalar(2)) / std::pow(_source_calibration[1], Scalar(2)) +
                      std::pow(_tmp7, Scalar(2)) / std::pow(_source_calibration[0], Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp10 = _tmp7 * _tmp9 / _source_calibration[0];
  const Scalar _tmp11 = -2 * std::pow(_source_pose[2], Scalar(2));
  const Scalar _tmp12 = -2 * std::pow(_source_pose[0], Scalar(2));
  const Scalar _tmp13 = _tmp8 * _tmp9 / _source_calibration[1];
  const Scalar _tmp14 = _source_pose[0] * _tmp5;
  const Scalar _tmp15 = _source_pose[2] * _tmp3;
  const Scalar _tmp16 = _tmp10 * (_tmp4 + _tmp6) + _tmp13 * (_tmp11 + _tmp12 + 1) +
                        _tmp9 * (-_tmp14 + _tmp15) +
                        source_inverse_range * (_source_pose[5] - _target_pose[5]);
  const Scalar _tmp17 = -2 * std::pow(_target_pose[2], Scalar(2));
  const Scalar _tmp18 = 1 - 2 * std::pow(_target_pose[1], Scalar(2));
  const Scalar _tmp19 = 1 - 2 * std::pow(_source_pose[1], Scalar(2));
  const Scalar _tmp20 = 2 * _source_pose[0] * _source_pose[2];
  const Scalar _tmp21 = _source_pose[1] * _tmp5;
  const Scalar _tmp22 = _tmp10 * (_tmp11 + _tmp19) + _tmp13 * (_tmp4 - _tmp6) +
                        _tmp9 * (_tmp20 + _tmp21) +
                        source_inverse_range * (_source_pose[4] - _target_pose[4]);
  const Scalar _tmp23 = _target_pose[0] * _tmp0;
  const Scalar _tmp24 = 2 * _target_pose[3];
  const Scalar _tmp25 = _target_pose[1] * _tmp24;
  const Scalar _tmp26 = _tmp10 * (_tmp20 - _tmp21) + _tmp13 * (_tmp14 + _tmp15) +
                        _tmp9 * (_tmp12 + _tmp19) +
                        source_inverse_range * (_source_pose[6] - _target_pose[6]);
  const Scalar _tmp27 = _target_pose[1] * _tmp0;
  const Scalar _tmp28 = _target_pose[0] * _tmp24;
  const Scalar _tmp29 = -2 * std::pow(_target_pose[0], Scalar(2));
  const Scalar _tmp30 =
      _tmp16 * (_tmp27 - _tmp28) + _tmp22 * (_tmp23 + _tmp25) + _tmp26 * (_tmp18 + _tmp29);
  const Scalar _tmp31 = Scalar(1.0) / (std::max<Scalar>(_tmp30, epsilon));

  // Output terms (2)
  if (reprojection_delta != nullptr) {
    Eigen::Matrix<Scalar, 2, 1>& _reprojection_delta = (*reprojection_delta);

    _reprojection_delta(0, 0) =
        _target_calibration[0] * _tmp31 *
            (_tmp16 * (_tmp1 + _tmp2) + _tmp22 * (_tmp17 + _tmp18) + _tmp26 * (_tmp23 - _tmp25)) +
        _target_calibration[2] - target_pixel(0, 0);
    _reprojection_delta(1, 0) = _target_calibration[1] * _tmp31 *
                                    (_tmp16 * (_tmp17 + _tmp29 + 1) + _tmp22 * (-_tmp1 + _tmp2) +
                                     _tmp26 * (_tmp27 + _tmp28)) +
                                _target_calibration[3] - target_pixel(1, 0);
  }

  if (is_valid != nullptr) {
    Scalar& _is_valid = (*is_valid);

    _is_valid = std::max<Scalar>(0, (((_tmp30) > 0) - ((_tmp30) < 0)));
  }
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym