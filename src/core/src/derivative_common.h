

#ifndef __DERIVATIVE_COMMON_H__
#define __DERIVATIVE_COMMON_H__

// Eigen
#include <Eigen/Core>

// Original
#include "rotations.h"

namespace core {

inline void ComputeInterimValues(const Eigen::Vector3d &T,
                                 const Eigen::Vector3d &rot,
                                 const Eigen::Matrix3d &K,
                                 const Eigen::Vector3d &p, Eigen::Matrix3d &R,
                                 Eigen::Matrix<double, 3, 4> &Proj,
                                 Eigen::Vector3d &x) {

  ConvertAngleAxisToRotationMatrix(rot, R);

  // X. Create camera.
  Proj.block<3, 3>(0, 0) = R;
  Proj.block<3, 1>(0, 3) = T;
  Proj = K * Proj;

  // X. Projection.
  x = Proj * p.homogeneous();
}

} // namespace core

#endif