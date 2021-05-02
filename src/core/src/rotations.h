// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2019 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: keir@google.com (Keir Mierle)
//         sameeragarwal@google.com (Sameer Agarwal)
//
// Templated functions for manipulating rotations. The templated
// functions are useful when implementing functors for automatic
// differentiation.
//
// In the following, the Quaternions are laid out as 4-vectors, thus:
//
//   q[0]  scalar part.
//   q[1]  coefficient of i.
//   q[2]  coefficient of j.
//   q[3]  coefficient of k.
//
// where: i*i = j*j = k*k = -1 and i*j = k, j*k = i, k*i = j.

#ifndef __ROTATIONS_H__
#define __ROTATIONS_H__

// STL
#include <numeric>

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace optimization {

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

} // namespace optimization

#endif