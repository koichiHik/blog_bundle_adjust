
// STL
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>

// Google
#include <gflags/gflags.h>
#include <glog/logging.h>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Original
#include "appl_utility.h"
#include <common_def.h>
#include <error_metric.h>
#include <optimizations.h>
#include <random_noise.h>
#include <utility.h>

DEFINE_string(camera_extrinsic_file_path, "", "");
DEFINE_string(camera_intrinsic_file_path, "", "");
DEFINE_string(image_location_file_path, "", "");
DEFINE_string(points_location_path, "", "");
DEFINE_string(result_save_dir, "", "");

using namespace core;

void ComputeCameraMatrix(const std::vector<Camera> &extrinsic_cams,
                         const std::map<size_t, size_t> &cam_intrinsic_map,
                         const std::vector<Eigen::Matrix3d> &intrinsic_cams,
                         std::vector<Camera> &cameras) {

  cameras.clear();
  for (size_t cam_idx = 0; cam_idx < extrinsic_cams.size(); cam_idx++) {
    size_t intrinsic_idx = cam_intrinsic_map.at(cam_idx);
    cameras.push_back(intrinsic_cams.at(intrinsic_idx) *
                      extrinsic_cams.at(cam_idx));
  }
}

int main(int argc, char **argv) {

  // X. Initial setting.
  google::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_alsologtostderr = 1;
  FLAGS_stderrthreshold = google::GLOG_INFO;
  google::InitGoogleLogging(argv[0]);

  // X. Load data.
  SfmData org_data;
  LoadData(FLAGS_camera_extrinsic_file_path, FLAGS_camera_intrinsic_file_path,
           FLAGS_image_location_file_path, FLAGS_points_location_path,
           org_data);

  // X. Noised data.
  SfmData noised_data;
  AddNoise(org_data, noised_data);

  // X. Compute Jacobian.
  std::vector<Camera> refined_extrinsics;
  std::vector<Eigen::Matrix3d> refined_intrinsics;
  std::vector<Eigen::Vector3d> refined_points3d;
  {
    core::CeresSolver().Optimize(noised_data.tracks, noised_data.extrinsic_cams,
                                 noised_data.intrinsic_cams,
                                 noised_data.extrinsic_intrinsic_map,
                                 noised_data.points3d, refined_extrinsics,
                                 refined_intrinsics, refined_points3d);
  }

  // X. Compute final reprojection error.
  {
    LOG(INFO) << "Reprojection Error (Final   Result) : "
              << core::ComputeReprojectionError(
                     noised_data.tracks, refined_intrinsics,
                     noised_data.extrinsic_intrinsic_map, refined_extrinsics,
                     refined_points3d);
    LOG(INFO) << "Average Reprojection Error (Final   Result) : "
              << core::ComputeAverageReprojectionError(
                     noised_data.tracks, refined_intrinsics,
                     noised_data.extrinsic_intrinsic_map, refined_extrinsics,
                     refined_points3d);
  }

  return 0;
}