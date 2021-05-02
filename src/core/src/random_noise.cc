
// X. STL
#include <random>

// X. Eigen
#include <Eigen/Geometry>

// X. Original
#include "random_noise.h"

namespace core {

void AddRandomNoiseToTrack(const std::vector<Track> &tracks, double mean,
                           double stddev, std::vector<Track> &noised_tracks) {

  std::default_random_engine engine(0);
  std::normal_distribution<> dist(mean, stddev);

  noised_tracks = tracks;
  for (size_t idx = 0; idx < noised_tracks.size(); idx++) {
    Track &noised_track = noised_tracks[idx];
    for (Track::iterator itr = noised_track.begin(); itr != noised_track.end();
         itr++) {
      itr->second(0) += dist(engine);
      itr->second(1) += dist(engine);
    }
  }
}

void AddRandomNoiseToPoints3d(const std::vector<Eigen::Vector3d> &points3d,
                              double mean, double stddev,
                              std::vector<Eigen::Vector3d> &noised_points3d) {

  std::default_random_engine engine(0);
  std::normal_distribution<> dist(mean, stddev);

  noised_points3d = points3d;
  for (size_t idx = 0; idx < noised_points3d.size(); idx++) {
    noised_points3d[idx](0) += dist(engine);
    noised_points3d[idx](1) += dist(engine);
    noised_points3d[idx](2) += dist(engine);
  }
}

void AddRandomNoiseToExtrinsic(const std::vector<Camera> &extrinsics,
                               double mean_trans, double stddev_trans,
                               double mean_angles, double stddev_angles,
                               std::vector<Camera> &noised_extrinsics) {

  std::default_random_engine engine(0);
  std::normal_distribution<> dists(mean_trans, stddev_trans);
  std::normal_distribution<> angles(mean_angles, stddev_angles);

  for (const Camera &extrinsic : extrinsics) {
    Eigen::Matrix3d R_noised;
    Eigen::Vector3d T_noised;

    Eigen::Vector3d euler = extrinsic.block<3, 3>(0, 0).eulerAngles(0, 1, 2);
    R_noised =
        Eigen::AngleAxisd(euler(0) + angles(engine), Eigen::Vector3d::UnitX()) *
        Eigen::AngleAxisd(euler(1) + angles(engine), Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(euler(2) + angles(engine), Eigen::Vector3d::UnitZ());

    T_noised(0) = extrinsic.block<3, 1>(0, 3)(0) + dists(engine);
    T_noised(1) = extrinsic.block<3, 1>(0, 3)(1) + dists(engine);
    T_noised(2) = extrinsic.block<3, 1>(0, 3)(2) + dists(engine);

    Camera noised_cam;
    noised_cam.block<3, 3>(0, 0) = R_noised;
    noised_cam.block<3, 1>(0, 3) = T_noised;
    noised_extrinsics.push_back(noised_cam);
  }
}

void AddRandomNoiseToIntrinsic(
    const std::vector<Eigen::Matrix3d> &intrinsics, double mean_fx,
    double stddev_fx, double mean_fy, double stddev_fy, double mean_cx,
    double stddev_cx, double mean_cy, double stddev_cy, double mean_skew,
    double stddev_skew, std::vector<Eigen::Matrix3d> &noised_intrinsics) {

  std::default_random_engine engine(0);
  std::normal_distribution<> dist_fx(mean_fx, stddev_fx);
  std::normal_distribution<> dist_fy(mean_fy, stddev_fy);
  std::normal_distribution<> dist_cx(mean_cx, stddev_cx);
  std::normal_distribution<> dist_cy(mean_cy, stddev_cy);
  std::normal_distribution<> dist_skew(mean_skew, stddev_skew);

  for (const Eigen::Matrix3d &K : intrinsics) {
    Eigen::Matrix3d K_ = K;
    K_(0, 0) += dist_fx(engine);
    K_(1, 1) += dist_fy(engine);
    K_(0, 2) += dist_cx(engine);
    K_(1, 2) += dist_cy(engine);
    K_(0, 1) += dist_skew(engine);
    noised_intrinsics.push_back(K_);
  }
}

} // namespace core