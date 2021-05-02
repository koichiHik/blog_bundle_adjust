
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
#include <file_utility.h>
#include <optimizations.h>

DEFINE_string(camera_extrinsic_file_path, "", "");
DEFINE_string(camera_intrinsic_file_path, "", "");
DEFINE_string(image_location_file_path, "", "");
DEFINE_string(points_location_path, "", "");
DEFINE_string(result_save_dir, "", "");

using namespace optimization;

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

  // X. Add noise.
  double image_scale, euclid_scale;
  SfmData noised_data;
  AddNoise(org_data, noised_data);

  // X. Normalize
  // NormalizeParameters(noised_data, image_scale, euclid_scale);

  // X. Compute Jacobian.
  SfmData refined_data(noised_data);
  {
    // optimization::GradientDescent().Optimize(
    optimization::GradientDescentWithLineSearch().Optimize(
        noised_data.tracks, noised_data.extrinsic_cams,
        noised_data.intrinsic_cams, noised_data.extrinsic_intrinsic_map,
        noised_data.points3d, refined_data.extrinsic_cams,
        refined_data.intrinsic_cams, refined_data.points3d);
  }

  // X. Re-Normalize
  // RenormalizeParameters(image_scale, euclid_scale, noised_data);
  // RenormalizeParameters(image_scale, euclid_scale, refined_data);

  // X. Save Data.
  SaveData(FLAGS_result_save_dir, "gd_noised", noised_data);
  SaveData(FLAGS_result_save_dir, "gd_refined", refined_data);

  // X. Compute final reprojection error.
  {
    LOG(INFO) << "Reprojection Error (Final   Result) : "
              << optimization::ComputeReprojectionError(
                     refined_data.tracks, refined_data.intrinsic_cams,
                     refined_data.extrinsic_intrinsic_map,
                     refined_data.extrinsic_cams, refined_data.points3d);
    LOG(INFO) << "Average Reprojection Error (Final   Result) : "
              << optimization::ComputeAverageReprojectionError(
                     refined_data.tracks, refined_data.intrinsic_cams,
                     refined_data.extrinsic_intrinsic_map,
                     refined_data.extrinsic_cams, refined_data.points3d);
  }

  return 0;
}