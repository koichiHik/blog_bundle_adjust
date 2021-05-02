

// STL
#include <vector>

// Boost
#include <boost/filesystem.hpp>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Google
#include <gflags/gflags.h>
#include <glog/logging.h>

// Original
#include <common_def.h>
#include <error_metric.h>
#include <optimizations.h>
#include <utility.h>

namespace optimization {

struct SfmData {
  std::vector<Camera> extrinsic_cams;
  std::map<size_t, size_t> extrinsic_intrinsic_map;
  std::vector<Eigen::Matrix3d> intrinsic_cams;
  std::vector<Eigen::Vector3d> points3d;
  std::vector<Track> tracks;
};

inline void AddRandomNoiseToTrack(const std::vector<Track> &tracks, double mean,
                                  double stddev,
                                  std::vector<Track> &noised_tracks) {

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

inline void
AddRandomNoiseToPoints3d(const std::vector<Eigen::Vector3d> &points3d,
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

inline void AddRandomNoiseToExtrinsic(const std::vector<Camera> &extrinsics,
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

inline void AddRandomNoiseToIntrinsic(
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

inline void LoadData(const std::string &extrinsic_filepath,
                     const std::string &intrinsic_filepath,
                     const std::string &image_location_filepath,
                     const std::string &points_location_filepath,
                     SfmData &org_data) {

  // X. Load extrinsics.
  utility::LoadCameraExtrinsicMatrix(extrinsic_filepath,
                                     org_data.extrinsic_cams);

  // X. Load intrinsics
  utility::LoadCameraIntrinsicMatrix(intrinsic_filepath,
                                     org_data.extrinsic_intrinsic_map,
                                     org_data.intrinsic_cams);

  // X. Load tracks.
  utility::LoadTracks(image_location_filepath, org_data.extrinsic_cams.size(),
                      org_data.tracks);

  // X. Load points 3d.
  utility::LoadPoints3D(points_location_filepath, org_data.points3d);

  // X. Create reprojection error of original and noised.
  LOG(INFO) << "Reprojection Error (Original Data) : "
            << optimization::ComputeReprojectionError(
                   org_data.tracks, org_data.intrinsic_cams,
                   org_data.extrinsic_intrinsic_map, org_data.extrinsic_cams,
                   org_data.points3d);
}

inline void SaveData(const std::string &result_directory,
                     const std::string &file_prefix, const SfmData &data) {

  boost::filesystem::path dir_path(result_directory);

  std::map<size_t, std::vector<size_t>> cam_indices_map;
  for (size_t idx = 0; idx < data.extrinsic_cams.size(); idx++) {
    const Camera &cam = data.extrinsic_cams[idx];
    utility::WriteCameraExtrinsicMatrix(
        dir_path.append(file_prefix + "_camera_extrinsics.txt")
            .generic_string(),
        cam.block<3, 3>(0, 0), cam.block<3, 1>(0, 3), false);
    cam_indices_map[data.extrinsic_intrinsic_map.at(idx)].push_back(idx);
  }

  for (size_t idx = 0; idx < data.intrinsic_cams.size(); idx++) {
    utility::WriteCameraIntrinsicMatrix(
        dir_path.append(file_prefix + "_camera_intrinsics.txt")
            .generic_string(),
        cam_indices_map[idx], data.intrinsic_cams[idx], false);
  }

  std::vector<std::vector<Eigen::Vector2d>> image_points;
  for (const Track &track : data.tracks) {

    std::vector<Eigen::Vector2d> points;
    for (Track::const_iterator citr = track.cbegin(); citr != track.cend();
         citr++) {
      points.push_back(citr->second);
    }
    image_points.push_back(points);
  }

  utility::WritePoints2D(
      dir_path.append(file_prefix + "_image_points2d.txt").generic_string(),
      image_points);

  utility::WritePoints3D(
      dir_path.append(file_prefix + "_refined_points3d.txt").generic_string(),
      data.points3d);
}

inline void AddNoise(const SfmData &org_data, SfmData &noised_data) {

  // X. Load extrinsics.
  // optimization::AddRandomNoiseToExtrinsic(org_data.extrinsic_cams, 0.0, 0.1,
  // 0.0,
  //                                0.05, noised_data.extrinsic_cams);
  optimization::AddRandomNoiseToExtrinsic(org_data.extrinsic_cams, 0.0, 0.0,
                                          0.0, 0.0, noised_data.extrinsic_cams);

  // X. Load intrinsics
  /*
  optimization::AddRandomNoiseToIntrinsic(org_data.intrinsic_cams, 0.0, 5.0,
  0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 0.0, noised_data.intrinsic_cams);
  */
  optimization::AddRandomNoiseToIntrinsic(org_data.intrinsic_cams, 0.0, 5.0,
                                          0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0,
                                          0.0, noised_data.intrinsic_cams);
  noised_data.extrinsic_intrinsic_map = org_data.extrinsic_intrinsic_map;

  // X. Load tracks.
  AddRandomNoiseToTrack(org_data.tracks, 0.0, 0.0, noised_data.tracks);

  // X. Load points 3d.
  // AddRandomNoiseToPoints3d(org_data.points3d, 0.0, 0.1,
  // noised_data.points3d);
  AddRandomNoiseToPoints3d(org_data.points3d, 0.0, 0.0, noised_data.points3d);

  // X. Create reprojection error of original and noised.
  LOG(INFO) << "Reprojection Error (Noised   Data) : "
            << optimization::ComputeReprojectionError(
                   noised_data.tracks, noised_data.intrinsic_cams,
                   noised_data.extrinsic_intrinsic_map,
                   noised_data.extrinsic_cams, noised_data.points3d);
}

inline void NormalizeParameters(SfmData &data, double &image_coord_scale,
                                double &euclid_coord_scale) {

  // X. Scaling for pixel unit.
  {
    double max_internal_p = 0;
    for (const Eigen::Matrix3d &K : data.intrinsic_cams) {
      max_internal_p = std::max(max_internal_p, K(0, 0));
      max_internal_p = std::max(max_internal_p, K(1, 1));
      max_internal_p = std::max(max_internal_p, K(0, 2));
      max_internal_p = std::max(max_internal_p, K(1, 2));
    }
    image_coord_scale = 1.0 / max_internal_p;
  }

  {
    for (Eigen::Matrix3d &K : data.intrinsic_cams) {
      K(0, 0) = K(0, 0) * image_coord_scale;
      K(1, 1) = K(1, 1) * image_coord_scale;
      K(0, 2) = K(0, 2) * image_coord_scale;
      K(1, 2) = K(1, 2) * image_coord_scale;
    }

    for (Track &track : data.tracks) {
      for (Track::iterator itr = track.begin(); itr != track.end(); itr++) {
        itr->second(0) = itr->second(0) * image_coord_scale;
        itr->second(1) = itr->second(1) * image_coord_scale;
      }
    }
  }

  // X. Scaling for euclid coordinate.
  {
    double max_coord = 0;
    for (const Eigen::Vector3d &p : data.points3d) {
      max_coord = std::max(std::abs(p(0)), max_coord);
      max_coord = std::max(std::abs(p(1)), max_coord);
      max_coord = std::max(std::abs(p(2)), max_coord);
    }

    for (const Camera cam : data.extrinsic_cams) {
      max_coord = std::max(std::abs(cam(0, 3)), max_coord);
      max_coord = std::max(std::abs(cam(1, 3)), max_coord);
      max_coord = std::max(std::abs(cam(2, 3)), max_coord);
    }
    euclid_coord_scale = 1.0 / max_coord;

    for (Eigen::Vector3d &p : data.points3d) {
      p(0) = p(0) * euclid_coord_scale;
      p(1) = p(1) * euclid_coord_scale;
      p(2) = p(2) * euclid_coord_scale;
    }

    for (Camera &cam : data.extrinsic_cams) {
      cam(0, 3) = cam(0, 3) * euclid_coord_scale;
      cam(1, 3) = cam(1, 3) * euclid_coord_scale;
      cam(2, 3) = cam(2, 3) * euclid_coord_scale;
    }
  }
}

inline void RenormalizeParameters(double image_scale, double euclid_scale,
                                  SfmData &data) {

  {
    for (Eigen::Matrix3d &K : data.intrinsic_cams) {
      K(0, 0) = K(0, 0) / image_scale;
      K(1, 1) = K(1, 1) / image_scale;
      K(0, 2) = K(0, 2) / image_scale;
      K(1, 2) = K(1, 2) / image_scale;
    }

    for (Track &track : data.tracks) {
      for (Track::iterator itr = track.begin(); itr != track.end(); itr++) {
        itr->second(0) = itr->second(0) / image_scale;
        itr->second(1) = itr->second(1) / image_scale;
      }
    }
  }

  for (Eigen::Vector3d &p : data.points3d) {
    p(0) = p(0) / euclid_scale;
    p(1) = p(1) / euclid_scale;
    p(2) = p(2) / euclid_scale;
  }

  for (Camera &cam : data.extrinsic_cams) {
    cam(0, 3) = cam(0, 3) / euclid_scale;
    cam(1, 3) = cam(1, 3) / euclid_scale;
    cam(2, 3) = cam(2, 3) / euclid_scale;
  }
}

} // namespace optimization