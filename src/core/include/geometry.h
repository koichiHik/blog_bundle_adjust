
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

// STL
#include <map>
#include <vector>

// Eigen
#include <Eigen/Core>

// Original
#include "common_def.h"

namespace core {

void DecomposeRQ(const Eigen::Matrix3d &M, Eigen::Matrix3d &R,
                 Eigen::Matrix3d &Q);

void DecomposeCameraMatrix(const Eigen::Matrix<double, 3, 4> &cam,
                           Eigen::Matrix3d &K, Eigen::Matrix3d &R,
                           Eigen::Vector3d &T);

void ComputeInternalCalibration(const Camera &cam, Eigen::Matrix3d &K,
                                Eigen::Matrix3d &R);

void ComputePoints3D(const std::vector<Track> &tracks,
                     const std::vector<Camera> &cams,
                     std::vector<Eigen::Vector3d> &tri_points);

Eigen::Matrix3d ComputeRotationMatrixViaEulerXYZ(double rx, double ry,
                                                 double rz);

std::vector<Eigen::Vector2d>
ProjectPoint(std::vector<Eigen::Vector3d> &points3d, Eigen::Matrix3d &R,
             Eigen::Vector3d &T, Eigen::Matrix3d &K);

} // namespace core

#endif