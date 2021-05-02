
#ifndef _VISUALIZE_H_
#define _VISUALIZE_H_

// STL
#include <string>
#include <vector>

// Eigen
#include <Eigen/Core>

// OpenCV
#include <opencv2/core.hpp>

// PCL
#include <pcl/PolygonMesh.h>

namespace utility {

void VisualizeWithPCL(const std::vector<Eigen::Vector3d> &points,
                      const std::vector<Eigen::Matrix<double, 3, 4>> &cams,
                      const std::vector<cv::Scalar> &colors);

void DrawPoints(const std::vector<Eigen::Vector2d> &image_points,
                const cv::Scalar &color, cv::Mat &image);

} // namespace utility

#endif