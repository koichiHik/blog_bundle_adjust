

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
#include <random_noise.h>
#include <utility.h>

namespace core {

struct SfmData {
  std::vector<Camera> extrinsic_cams;
  std::map<size_t, size_t> extrinsic_intrinsic_map;
  std::vector<Eigen::Matrix3d> intrinsic_cams;
  std::vector<Eigen::Vector3d> points3d;
  std::vector<Track> tracks;
};

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
            << core::ComputeReprojectionError(
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
  // core::AddRandomNoiseToExtrinsic(org_data.extrinsic_cams, 0.0, 0.1, 0.0,
  //                                0.05, noised_data.extrinsic_cams);
  core::AddRandomNoiseToExtrinsic(org_data.extrinsic_cams, 0.0, 0.0, 0.0, 0.0,
                                  noised_data.extrinsic_cams);

  // X. Load intrinsics
  /*
  core::AddRandomNoiseToIntrinsic(org_data.intrinsic_cams, 0.0, 5.0, 0.0, 5.0,
                                  0.0, 5.0, 0.0, 5.0, 0.0, 0.0,
                                  noised_data.intrinsic_cams);
  */
  core::AddRandomNoiseToIntrinsic(org_data.intrinsic_cams, 0.0, 5.0, 0.0, 5.0,
                                  0.0, 5.0, 0.0, 5.0, 0.0, 0.0,
                                  noised_data.intrinsic_cams);
  noised_data.extrinsic_intrinsic_map = org_data.extrinsic_intrinsic_map;

  // X. Load tracks.
  AddRandomNoiseToTrack(org_data.tracks, 0.0, 0.0, noised_data.tracks);

  // X. Load points 3d.
  // AddRandomNoiseToPoints3d(org_data.points3d, 0.0, 0.1,
  // noised_data.points3d);
  AddRandomNoiseToPoints3d(org_data.points3d, 0.0, 0.0, noised_data.points3d);

  // X. Create reprojection error of original and noised.
  LOG(INFO) << "Reprojection Error (Noised   Data) : "
            << core::ComputeReprojectionError(
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

} // namespace core