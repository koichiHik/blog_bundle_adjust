
#ifndef __JACOBI_COMMON_H__
#define __JACOBI_COMMON_H__

// STL
#include <numeric>

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace core {

enum HandedNess {
  Right,
  Left,
};

inline void RotationMatrixToQuaternion(const Eigen::Matrix3d &R,
                                       Eigen::Vector4d &quaternion) {
  const double trace = R(0, 0) + R(1, 1) + R(2, 2);
  if (trace >= 0.0) {
    double t = sqrt(trace + 1.0);
    quaternion(0) = 0.5 * t;
    t = 0.5 / t;
    quaternion(1) = (R(2, 1) - R(1, 2)) * t;
    quaternion(2) = (R(0, 2) - R(2, 0)) * t;
    quaternion(3) = (R(1, 0) - R(0, 1)) * t;
  } else {
    int i = 0;
    if (R(1, 1) > R(0, 0)) {
      i = 1;
    }

    if (R(2, 2) > R(i, i)) {
      i = 2;
    }

    const int j = (i + 1) % 3;
    const int k = (j + 1) % 3;
    double t = sqrt(R(i, i) - R(j, j) - R(k, k) + 1.0);
    quaternion(i + 1) = 0.5 * t;
    t = 0.5 / t;
    quaternion(0) = (R(k, j) - R(j, k)) * t;
    quaternion(j + 1) = (R(j, i) + R(i, j)) * t;
    quaternion(k + 1) = (R(k, i) + R(i, k)) * t;
  }
}

inline void QuaternionToAngleAxis(const Eigen::Vector4d &quaternion,
                                  Eigen::Vector3d &angle_axis) {
  const double &q1 = quaternion(1);
  const double &q2 = quaternion(2);
  const double &q3 = quaternion(3);
  const double sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;

  // For quaternions representing non-zero rotation, the conversion
  // is numerically stable.
  if (sin_squared_theta > 0.0) {
    const double sin_theta = sqrt(sin_squared_theta);
    const double &cos_theta = quaternion(0);

    // If cos_theta is negative, theta is greater than pi/2, which
    // means that angle for the angle_axis vector which is 2 * theta
    // would be greater than pi.
    //
    // While this will result in the correct rotation, it does not
    // result in a normalized angle-axis vector.
    //
    // In that case we observe that 2 * theta ~ 2 * theta - 2 * pi,
    // which is equivalent saying
    //
    //   theta - pi = atan(sin(theta - pi), cos(theta - pi))
    //              = atan(-sin(theta), -cos(theta))
    //
    const double two_theta =
        2.0 * ((cos_theta < 0.0) ? atan2(-sin_theta, -cos_theta)
                                 : atan2(sin_theta, cos_theta));
    const double k = two_theta / sin_theta;
    angle_axis(0) = q1 * k;
    angle_axis(1) = q2 * k;
    angle_axis(2) = q3 * k;
  } else {
    // For zero rotation, sqrt() will produce NaN in the derivative since
    // the argument is zero.  By approximating with a Taylor series,
    // and truncating at one term, the value and first derivatives will be
    // computed correctly when Jets are used.
    const double k = 2.0;
    angle_axis(0) = q1 * k;
    angle_axis(1) = q2 * k;
    angle_axis(2) = q3 * k;
  }
}

inline void ConvertRotationMatrixToAngleAxis(const Eigen::Matrix3d &R,
                                             Eigen::Vector3d &angle_axis) {

  Eigen::Vector4d quaternion;
  RotationMatrixToQuaternion(R, quaternion);
  QuaternionToAngleAxis(quaternion, angle_axis);
}

inline void ConvertAngleAxisToRotationMatrix(const Eigen::Vector3d &angle_axis,
                                             Eigen::Matrix3d &R) {

  static const double kOne = 1.0;

  double theta2 = angle_axis.dot(angle_axis);
  if (theta2 > std::numeric_limits<double>::epsilon()) {
    double theta = std::sqrt(theta2);
    double wx = angle_axis(0) / theta;
    double wy = angle_axis(1) / theta;
    double wz = angle_axis(2) / theta;

    double costheta = std::cos(theta);
    double sintheta = std::sin(theta);

    R(0, 0) = costheta + wx * wx * (kOne - costheta);
    R(1, 0) = wz * sintheta + wx * wy * (kOne - costheta);
    R(2, 0) = -wy * sintheta + wx * wz * (kOne - costheta);
    R(0, 1) = wx * wy * (kOne - costheta) - wz * sintheta;
    R(1, 1) = costheta + wy * wy * (kOne - costheta);
    R(2, 1) = wx * sintheta + wy * wz * (kOne - costheta);
    R(0, 2) = wy * sintheta + wx * wz * (kOne - costheta);
    R(1, 2) = -wx * sintheta + wy * wz * (kOne - costheta);
    R(2, 2) = costheta + wz * wz * (kOne - costheta);

  } else {

    R(0, 0) = kOne;
    R(1, 0) = angle_axis(2);
    R(2, 0) = -angle_axis(1);
    R(0, 1) = -angle_axis(2);
    R(1, 1) = kOne;
    R(2, 1) = angle_axis(0);
    R(0, 2) = angle_axis(1);
    R(1, 2) = -angle_axis(0);
    R(2, 2) = kOne;
  }
}

/*
inline void ConvertAngleAxisToRotationMatrix(const Eigen::Vector3d &angle_axis,
                                             Eigen::Matrix3d &R) {

  double theta2 = angle_axis.dot(angle_axis);
  Eigen::AngleAxisd eigen_angle_axis(std::sqrt(theta2),
                                     angle_axis.normalized());
  R = eigen_angle_axis.toRotationMatrix();
}
*/

} // namespace core

#endif