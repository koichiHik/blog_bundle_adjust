

// STL
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Google
#include <gflags/gflags.h>
#include <glog/logging.h>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// OpenCV
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// Original
#include <file_utility.h>
#include <geometry.h>
#include <visualize.h>

DEFINE_string(camera_extrinsic_file_path, "", "");
DEFINE_string(camera_intrinsic_file_path, "", "");
DEFINE_string(image_location_file_path, "", "");
DEFINE_string(points_location_path, "", "");

std::vector<Eigen::Vector3d> CreateArtificialPoint3D() {

  size_t num = 30.0;
  double R = 2;
  double height = 1;
  double dist = 4;
  double d_theta = 2 * M_PI / static_cast<double>(num);
  double theta = M_PI;

  std::vector<Eigen::Vector3d> points3d;
  for (size_t idx = 0; idx < num; idx++) {
    double x = R / 2.0 * cos(d_theta * idx);
    double y = R / 2.0 * sin(d_theta * idx) - height;
    double z = dist;
    Eigen::Vector3d p(x, y, z);
    points3d.push_back(p);
  }

  return points3d;
}

int main(int argc, char **argv) {

  // X. Initial setting.
  google::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_alsologtostderr = 1;
  FLAGS_stderrthreshold = google::GLOG_INFO;
  google::InitGoogleLogging(argv[0]);

  std::vector<Eigen::Vector3d> points3d = CreateArtificialPoint3D();

  // X. Create intrinsics.
  Eigen::Matrix3d K1, K2;
  K1 << 700, 0, 500, 0, 700, 500, 0, 0, 1;
  K2 << 350, 0, 500, 0, 350, 500, 0, 0, 1;

  // X. Create extrinsics.
  Eigen::Matrix3d R1_cam_to_org, R2_cam_to_org, R3_cam_to_org;
  Eigen::Vector3d T1_cam_to_org, T2_cam_to_org, T3_cam_to_org;
  std::vector<Eigen::Vector2d> image_points2d_1, image_points2d_2,
      image_points2d_3;
  {
    double tx = 0, ty = 0, tz = 0;
    double rx = 0, ry = 0, rz = 0;
    R1_cam_to_org = optimization::ComputeRotationMatrixViaEulerXYZ(
                        rx / 180.0 * M_PI, ry / 180.0 * M_PI, rz / 180.0 * M_PI)
                        .transpose();
    T1_cam_to_org = -R1_cam_to_org * Eigen::Vector3d(tx, ty, tz);
    image_points2d_1 =
        optimization::ProjectPoint(points3d, R1_cam_to_org, T1_cam_to_org, K1);
  }

  {
    double tx = 3, ty = 0, tz = 0;
    double rx = 0, ry = -30, rz = 0;
    R2_cam_to_org = optimization::ComputeRotationMatrixViaEulerXYZ(
                        rx / 180.0 * M_PI, ry / 180.0 * M_PI, rz / 180.0 * M_PI)
                        .transpose();
    T2_cam_to_org = -R2_cam_to_org * Eigen::Vector3d(tx, ty, tz);
    image_points2d_2 =
        optimization::ProjectPoint(points3d, R2_cam_to_org, T2_cam_to_org, K1);
  }

  {
    double tx = -3.0, ty = 0, tz = 1.0;
    double rx = 0, ry = 15, rz = 0;
    R3_cam_to_org = optimization::ComputeRotationMatrixViaEulerXYZ(
                        rx / 180.0 * M_PI, ry / 180.0 * M_PI, rz / 180.0 * M_PI)
                        .transpose();
    T3_cam_to_org = -R3_cam_to_org * Eigen::Vector3d(tx, ty, tz);
    image_points2d_3 =
        optimization::ProjectPoint(points3d, R3_cam_to_org, T3_cam_to_org, K2);
  }

  std::vector<cv::Scalar> colors{cv::Scalar(0, 0, 255), cv::Scalar(0, 255, 0),
                                 cv::Scalar(255, 0, 0)};
  cv::Mat img1(1000, 1000, CV_8UC3), img2(1000, 1000, CV_8UC3),
      img3(1000, 1000, CV_8UC3);
  img1 = cv::Scalar::all(255);
  img2 = cv::Scalar::all(255);
  img3 = cv::Scalar::all(255);
  utility::DrawPoints(image_points2d_1, colors[0], img1);
  utility::DrawPoints(image_points2d_2, colors[1], img2);
  utility::DrawPoints(image_points2d_3, colors[2], img3);

  cv::imshow("Image of Camera 1", img1);
  cv::imshow("Image of Camera 2", img2);
  cv::imshow("Image of Camera 3", img3);
  cv::waitKey(100);

  {
    std::vector<Eigen::Matrix<double, 3, 4>> cams(3);
    cams[0].block<3, 3>(0, 0) = R1_cam_to_org;
    cams[0].block<3, 1>(0, 3) = T1_cam_to_org;
    cams[1].block<3, 3>(0, 0) = R2_cam_to_org;
    cams[1].block<3, 1>(0, 3) = T2_cam_to_org;
    cams[2].block<3, 3>(0, 0) = R3_cam_to_org;
    cams[2].block<3, 1>(0, 3) = T3_cam_to_org;
    utility::VisualizeWithPCL(points3d, cams, colors);
  }

  // X. Save data.
  {
    // X. Save image points.
    utility::WritePoints2D(
        FLAGS_image_location_file_path,
        std::vector<std::vector<Eigen::Vector2d>>{
            image_points2d_1, image_points2d_2, image_points2d_3});

    // X. Save camera extrinsics.
    utility::WriteCameraExtrinsicMatrix(FLAGS_camera_extrinsic_file_path,
                                        R1_cam_to_org, T1_cam_to_org, true);
    utility::WriteCameraExtrinsicMatrix(FLAGS_camera_extrinsic_file_path,
                                        R2_cam_to_org, T2_cam_to_org, false);
    utility::WriteCameraExtrinsicMatrix(FLAGS_camera_extrinsic_file_path,
                                        R3_cam_to_org, T3_cam_to_org, false);

    // X. Save camera intrinsics.
    utility::WriteCameraIntrinsicMatrix(FLAGS_camera_intrinsic_file_path,
                                        std::vector<size_t>{0, 1}, K1, true);
    utility::WriteCameraIntrinsicMatrix(FLAGS_camera_intrinsic_file_path,
                                        std::vector<size_t>{2}, K2, false);

    // X. Save point location.
    utility::WritePoints3D(FLAGS_points_location_path, points3d);
  }

  return 0;
}