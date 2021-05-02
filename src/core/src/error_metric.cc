

// X. Eigen
#include <Eigen/Dense>

// X. Google
#include <glog/logging.h>

// X. Original
#include "differential_operators.h"
#include "rotations.h"
#include <error_metric.h>

namespace core {

double ComputeReprojectionErrorWithStepLength(
    const Eigen::MatrixXd &gradient, double step_length,
    const std::vector<Track> &tracks_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map_src) {

  // X. Copy
  std::vector<Eigen::Vector3d> T_dst = T_src;
  std::vector<Eigen::Matrix3d> K_dst = K_src;
  std::vector<Eigen::Vector3d> Rot_dst = Rot_src;
  std::vector<Eigen::Vector3d> points3d_dst = points3d_src;

  // X. Update parameters.
  UpdateParameters(step_length * gradient, tracks_src,
                   extrinsic_intrinsic_map_src, K_dst, T_dst, Rot_dst,
                   points3d_dst);

  return ComputeReprojectionError(tracks_src, K_dst,
                                  extrinsic_intrinsic_map_src, Rot_dst, T_dst,
                                  points3d_dst);
}

double ComputeReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Eigen::Vector3d> &Rot,
    const std::vector<Eigen::Vector3d> &T,
    const std::vector<Eigen::Vector3d> &points3d) {

  std::vector<Camera> extrinsics;
  for (size_t cam_ext_idx = 0; cam_ext_idx < T.size(); cam_ext_idx++) {
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);

    Camera extrinsic;
    extrinsic.block<3, 3>(0, 0) = R;
    extrinsic.block<3, 1>(0, 3) = T[cam_ext_idx];
    extrinsics.push_back(extrinsic);
  }

  return ComputeReprojectionError(tracks, intrinsics, extrinsic_intrinsic_map,
                                  extrinsics, points3d);
}

double ComputeReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Vector3d> &points3d) {

  double reproj_err = 0;
  for (size_t idx = 0; idx < tracks.size(); idx++) {
    const Track &track = tracks[idx];
    for (Track::const_iterator citr = track.cbegin(); citr != track.cend();
         citr++) {
      reproj_err +=
          0.5 *
          std::pow(((*citr).second -
                    (intrinsics[extrinsic_intrinsic_map.at(citr->first)] *
                     extrinsics[citr->first] * points3d[idx].homogeneous())
                        .hnormalized())
                       .norm(),
                   2.0);
    }
  }
  return reproj_err;
}

double ComputeAverageReprojectionError(
    const std::vector<Track> &tracks,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Vector3d> &points3d) {
  double reproj_err = ComputeReprojectionError(
      tracks, intrinsics, extrinsic_intrinsic_map, extrinsics, points3d);

  size_t num = 0;
  for (size_t idx = 0; idx < tracks.size(); idx++) {
    const Track &track = tracks[idx];
    for (Track::const_iterator citr = track.cbegin(); citr != track.cend();
         citr++) {
      num++;
    }
  }
  return reproj_err / num;
}

} // namespace core