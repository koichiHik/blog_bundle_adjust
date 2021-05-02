

// Google
#include <glog/logging.h>

// OpenCV
#include <opencv2/imgproc.hpp>

// PCL
#include <pcl/ros/conversions.h>
#include <pcl/visualization/pcl_visualizer.h>

// Original
#include "geometry.h"
#include "visualize.h"

namespace utility {

float CAM_SCALE_COEFF = 0.25;
int CAM_POLYGON[18] = {0, 1, 2, 0, 3, 1, 0, 4, 3, 0, 2, 4, 3, 1, 4, 2, 4, 1};

pcl::PointXYZRGB EigenToPointXYZRGB(Eigen::Vector3d v, Eigen::Vector3d rgb) {
  pcl::PointXYZRGB p(rgb[0], rgb[1], rgb[2]);
  p.x = v[0];
  p.y = v[1];
  p.z = v[2];
  return p;
}

void ComposeCameraElement(
    const Eigen::Matrix3d &R, const Eigen::Vector3d &T, const std::string &name,
    float r, float g, float b,
    std::vector<std::pair<std::string, pcl::PolygonMesh>> &cam_meshes,
    double s) {

  Eigen::Vector3d euler = R.eulerAngles(0, 1, 2);
  Eigen::Vector3d T_trans = -R.transpose() * T;
  Eigen::Vector3d v_right, v_up, v_forward;
  v_right = R.row(0).normalized() * s;
  v_up = -R.row(1).normalized() * s;
  v_forward = R.row(2).normalized() * s;
  Eigen::Vector3d rgb(r, g, b);

  // Polygon Mesh
  {
    pcl::PointCloud<pcl::PointXYZRGB> mesh_cld;
    mesh_cld.push_back(EigenToPointXYZRGB(T_trans, rgb));
    mesh_cld.push_back(EigenToPointXYZRGB(
        T_trans + v_forward + v_right / 2.0 + v_up / 2.0, rgb));
    mesh_cld.push_back(EigenToPointXYZRGB(
        T_trans + v_forward + v_right / 2.0 - v_up / 2.0, rgb));
    mesh_cld.push_back(EigenToPointXYZRGB(
        T_trans + v_forward - v_right / 2.0 + v_up / 2.0, rgb));
    mesh_cld.push_back(EigenToPointXYZRGB(
        T_trans + v_forward - v_right / 2.0 - v_up / 2.0, rgb));

    pcl::PolygonMesh pm;
    pm.polygons.resize(6);
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++) {
        pm.polygons[i].vertices.push_back(CAM_POLYGON[i * 3 + j]);
      }
    }
    pcl::toROSMsg(mesh_cld, pm.cloud);
    cam_meshes.push_back(std::make_pair(name, pm));
  }
}

void ConvertFromOpenCVToPCL(const Eigen::Matrix3d &R_cv, Eigen::Matrix3d &R_pcl,
                            const Eigen::Vector3d &T_cv,
                            Eigen::Vector3d &T_pcl) {

  Eigen::Vector3d euler_cv = R_cv.eulerAngles(0, 1, 2);
  R_pcl = optimization::ComputeRotationMatrixViaEulerXYZ(
      -euler_cv(0), -euler_cv(1), euler_cv(2));
  T_pcl(0) = -T_cv(0);
  T_pcl(1) = -T_cv(1);
  T_pcl(2) = T_cv(2);
}

void VisualizeWithPCL(const std::vector<Eigen::Vector3d> &points,
                      const std::vector<Eigen::Matrix<double, 3, 4>> &cams,
                      const std::vector<cv::Scalar> &colors) {

  // X. Visualization via PCL.
  pcl::visualization::PCLVisualizer viewer("pcl_viewer");
  viewer.addCoordinateSystem(0.5, 0);
  viewer.initCameraParameters();
  viewer.setCameraPosition(0, 0, -7.5, 0, 0, 0, 0, 0, 0);

  // X. Convert point.
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud(
      new pcl::PointCloud<pcl::PointXYZRGB>());
  for (const Eigen::Vector3d &pnt_cv : points) {

    pcl::PointXYZRGB pcl_p;
    pcl_p.x = -pnt_cv(0);
    pcl_p.y = -pnt_cv(1);
    pcl_p.z = pnt_cv(2);
    pcl_p.r = 200;
    pcl_p.g = 200;
    pcl_p.b = 0;

    point_cloud->push_back(pcl_p);
  }

  // X. Add point cloud.
  viewer.addPointCloud(point_cloud, "cloud");
  viewer.setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "cloud");

  // X. Add camera frustum.
  for (size_t idx = 0; idx < cams.size(); idx++) {
    std::string name = "cam" + std::to_string(idx);
    std::vector<std::pair<std::string, pcl::PolygonMesh>> cam_meshes;

    Eigen::Matrix3d R_cv = cams[idx].block<3, 3>(0, 0);
    Eigen::Vector3d T_cv = cams[idx].block<3, 1>(0, 3);
    cv::Scalar color = colors[idx];

    Eigen::Matrix3d R_pcl;
    Eigen::Vector3d T_pcl;
    ConvertFromOpenCVToPCL(R_cv, R_pcl, T_cv, T_pcl);
    ComposeCameraElement(R_pcl, T_pcl, name, color(2), color(1), color(0),
                         cam_meshes, CAM_SCALE_COEFF);
    viewer.addPolygonMesh(cam_meshes.back().second, cam_meshes.back().first);
  }

  while (!viewer.wasStopped()) {
    viewer.spinOnce();
  }
}

void DrawPoints(const std::vector<Eigen::Vector2d> &image_points,
                const cv::Scalar &color, cv::Mat &image) {

  for (size_t idx = 0; idx < image_points.size(); idx++) {
    int x = image_points[idx](0);
    int y = image_points[idx](1);
    cv::circle(image, cv::Point(x, y), 2, color, 2);
  }
}

} // namespace utility
