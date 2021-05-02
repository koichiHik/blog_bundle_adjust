
// STL
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <thread>

// Google
#include <gflags/gflags.h>
#include <glog/logging.h>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

// OpenCV
#include <opencv2/calib3d.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>

// PCL
#include <pcl/visualization/pcl_visualizer.h>

// Original
#include <bundle_adjust.h>
#include <common_def.h>
#include <geometry.h>
#include <utility.h>

DEFINE_string(camera_file_path, "", "");
DEFINE_string(image_location_file_path, "", "");
DEFINE_string(points_location_path, "", "");

using namespace core;

int main(int argc, char **argv) {
  // X. Initial setting.
  google::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_alsologtostderr = 1;
  FLAGS_stderrthreshold = google::GLOG_INFO;
  google::InitGoogleLogging(argv[0]);

  LOG(INFO) << "Compute track locations";
  LOG(INFO) << "Camera File Path : " << FLAGS_camera_file_path;
  LOG(INFO) << "Image Location File Path : " << FLAGS_image_location_file_path;

#if 0

  std::vector<Camera> cameras;
  utility::LoadCameraMatrixMap(FLAGS_camera_file_path, cameras);

  std::vector<Camera> mod_cameras = cameras;
  for (int idx = 0; idx < cameras.size(); idx++) {
    Eigen::Matrix3d K, R;
    Eigen::Vector3d T;
    DecomposeCameraMatrix(cameras[idx], K, R, T);
    mod_cameras[idx].block<3, 3>(0, 0) = K * R;
    mod_cameras[idx].block<3, 1>(0, 3) = K * T;
  }
  cameras = mod_cameras;

  std::vector<Track> tracks;
  utility::LoadTracks(FLAGS_image_location_file_path, cameras.size(), tracks);

  std::vector<Eigen::Vector3d> tri_points;
  ComputePoints3D(tracks, cameras, tri_points);

  std::vector<Camera> refined_cameras;
  std::vector<Eigen::Vector3d> refined_points;
  core::BundleAdjuster(cameras, tracks, tri_points, refined_cameras,
                       refined_points);

  utility::WritePoints3D(FLAGS_points_location_path, tri_points);

  pcl::PointCloud<pcl::PointXYZRGB> cloud;
  for (const auto &p : tri_points) {

    pcl::PointXYZRGB pcl_p;
    pcl_p.x = p(0);
    pcl_p.y = p(1);
    pcl_p.z = p(2);
    pcl_p.r = 250;
    pcl_p.g = 250;
    pcl_p.b = 250;
    cloud.push_back(pcl_p);
  }

  pcl::visualization::PCLVisualizer viewer("Dinosaur");
  viewer.addPointCloud(cloud.makeShared(), "cloud");
  viewer.addCoordinateSystem(0.2);
  viewer.spin();
#endif

  return 0;
}
