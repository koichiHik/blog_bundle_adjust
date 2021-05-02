
#ifndef __RANDOM_NOISE_H__
#define __RANDOM_NOISE_H__

// X. STL
#include <vector>

// X. Original
#include <common_def.h>

namespace core {

void AddRandomNoiseToTrack(const std::vector<Track> &tracks, double mean,
                           double stddev, std::vector<Track> &noised_tracks);

void AddRandomNoiseToPoints3d(const std::vector<Eigen::Vector3d> &points3d,
                              double mean, double stddev,
                              std::vector<Eigen::Vector3d> &noised_points3d);

void AddRandomNoiseToExtrinsic(const std::vector<Camera> &extrinsics,
                               double mean_trans, double stddev_trans,
                               double mean_angles, double stddev_angles,
                               std::vector<Camera> &noised_extrinsics);

void AddRandomNoiseToIntrinsic(const std::vector<Eigen::Matrix3d> &intrinsics,
                               double mean_fx, double stddev_fx, double mean_fy,
                               double stddev_fy, double mean_cx,
                               double stddev_cx, double mean_cy,
                               double stddev_cy, double mean_skew,
                               double stddev_skew,
                               std::vector<Eigen::Matrix3d> &noised_intrinsics);

} // namespace core

#endif
