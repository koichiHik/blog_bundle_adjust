
#ifndef __FILE_UTIL_H__
#define __FILE_UTIL_H__

// STL
#include <string>
#include <vector>

// Eigen
#include <Eigen/Core>

// X. Original
#include <common_def.h>

namespace utility {

void WritePoints3D(const std::string &path,
                   const std::vector<Eigen::Vector3d> &points3d);

void WritePoints2D(const std::string &path,
                   const std::vector<std::vector<Eigen::Vector2d>> &points2ds);

void WriteCameraExtrinsicMatrix(const std::string &path,
                                const Eigen::Matrix3d &R,
                                const Eigen::Vector3d &T, bool overwrite);

void WriteCameraIntrinsicMatrix(const std::string &path,
                                std::vector<size_t> cam_indices,
                                const Eigen::Matrix3d &K, bool overwrite);

void LoadPoints3D(const std::string &path,
                  std::vector<Eigen::Vector3d> &points3d);

void LoadCameraIntrinsicMatrix(const std::string &path,
                               std::map<size_t, size_t> &cam_intrinsic_map,
                               std::vector<Eigen::Matrix3d> &cam_intrinsics);

void LoadCameraExtrinsicMatrix(
    const std::string &path, std::vector<Eigen::Matrix<double, 3, 4>> &cameras);

void LoadTracks(const std::string &path, int cam_no,
                std::vector<optimization::Track> &tracks);

} // namespace utility

#endif