
#ifndef __ERROR_METRIC_H__
#define __ERROR_METRIC_H__

// X. STL
#include <vector>

// X. Eigen
#include <Eigen/Core>

// X. Original
#include <common_def.h>

namespace optimization {

double ComputeReprojectionErrorWithStepLength(
    const Eigen::MatrixXd &gradient, double step_length,
    const std::vector<Track> &tracks_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map_src);

double ComputeReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Eigen::Vector3d> &R,
    const std::vector<Eigen::Vector3d> &T,
    const std::vector<Eigen::Vector3d> &points3d);

double ComputeReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Vector3d> &points3d);

double ComputeAverageReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Vector3d> &points3d);

} // namespace optimization

#endif