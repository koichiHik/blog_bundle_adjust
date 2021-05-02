
// X. STL
#include <iostream>

// X. Google
#include <glog/logging.h>

// X. Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

// X. OpenCV
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc.hpp>

// X. Original
#include "geometry.h"

namespace optimization {

void PrintMatrix(const Eigen::Matrix3d &m) {

  for (int row = 0; row < m.cols(); row++) {
    for (int col = 0; col < m.rows(); col++) {
      std::cout << m(row, col) << ",";
    }
    std::cout << std::endl;
  }
}

void DecomposeRQ(const Eigen::Matrix3d &M, Eigen::Matrix3d &R,
                 Eigen::Matrix3d &Q) {

  R = M;
  Eigen::Matrix3d Rx, Ry, Rz;

  // X. Compute X.
  {
    double sx = R(2, 1) / std::sqrt(R(2, 2) * R(2, 2) + R(2, 1) * R(2, 1));
    double cx = -R(2, 2) / std::sqrt(R(2, 2) * R(2, 2) + R(2, 1) * R(2, 1));
    Rx << 1, 0, 0, 0, cx, -sx, 0, sx, cx;
    R = R * Rx;
  }

  // X. Compute Y.
  {
    double sy = R(2, 0) / std::sqrt(R(2, 2) * R(2, 2) + R(2, 0) * R(2, 0));
    double cy = R(2, 2) / std::sqrt(R(2, 2) * R(2, 2) + R(2, 0) * R(2, 0));
    Ry << cy, 0, sy, 0, 1, 0, -sy, 0, cy;
    R = R * Ry;
  }

  // X. Compute Z.
  {
    double sz = R(1, 0) / std::sqrt(R(1, 1) * R(1, 1) + R(1, 0) * R(1, 0));
    double cz = -R(1, 1) / std::sqrt(R(1, 1) * R(1, 1) + R(1, 0) * R(1, 0));
    Rz << cz, -sz, 0, sz, cz, 0, 0, 0, 1;
    R = R * Rz;
  }

  Q = Rz.transpose() * Ry.transpose() * Rx.transpose();

  // X. Force diagonal entry of R positive.
  Eigen::Vector3d sign = R.diagonal().cwiseSign();
  Eigen::DiagonalMatrix<double, 3> diag(sign);
  R = R * diag;
  Q = diag * Q;
}

void DecomposeCameraMatrix(const Eigen::Matrix<double, 3, 4> &cam,
                           Eigen::Matrix3d &K, Eigen::Matrix3d &R,
                           Eigen::Vector3d &T) {

  Eigen::Matrix3d M = cam.block<3, 3>(0, 0);
  DecomposeRQ(M, K, R);
  T = K.inverse() * cam.block<3, 1>(0, 3);
}

void ComputeInternalCalibration(const Eigen::Matrix<double, 3, 4> &cam,
                                Eigen::Matrix3d &K, Eigen::Matrix3d &R) {

  Eigen::Matrix3d M = cam.block<3, 3>(0, 0);
  DecomposeRQ(M, K, R);
}

void CreatePairSet(std::vector<int> &single_keys,
                   std::vector<std::pair<int, int>> &pair_keys) {

  std::sort(single_keys.begin(), single_keys.end());

  pair_keys.clear();
  for (size_t idx1 = 0; idx1 < single_keys.size() - 1; idx1++) {
    for (size_t idx2 = idx1 + 1; idx2 < single_keys.size(); idx2++) {
      pair_keys.push_back(std::make_pair(single_keys[idx1], single_keys[idx2]));
    }
  }
}

Eigen::Vector3d TriangulateTrack(const Track &track,
                                 const std::vector<Camera> &cams) {

  // X. Create pair set that sees this track.
  std::vector<std::pair<int, int>> pair_keys;
  {
    std::vector<int> single_keys;
    for (const auto e : track) {
      single_keys.push_back(e.first);
    }
    CreatePairSet(single_keys, pair_keys);
  }

  // X. Create linear matrix for triangulation.
  Eigen::Matrix<double, Eigen::Dynamic, 4> A(4 * pair_keys.size(), 4);
  for (size_t idx = 0; idx < pair_keys.size(); idx++) {

    CHECK(pair_keys[idx].first < static_cast<int>(cams.size()));
    CHECK(pair_keys[idx].second < static_cast<int>(cams.size()));
    CHECK(track.count(pair_keys[idx].first));
    CHECK(track.count(pair_keys[idx].second));

    const auto &cam1 = cams[pair_keys[idx].first];
    const auto &p1 = track.at(pair_keys[idx].first);
    const auto &cam2 = cams[pair_keys[idx].second];
    const auto &p2 = track.at(pair_keys[idx].second);

    A.row(4 * idx) = p1(0) * cam1.row(2) - cam1.row(0);
    A.row(4 * idx + 1) = p1(1) * cam1.row(2) - cam1.row(1);
    A.row(4 * idx + 2) = p2(0) * cam2.row(2) - cam2.row(0);
    A.row(4 * idx + 3) = p2(1) * cam2.row(2) - cam2.row(1);
  }

  Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> svd(
      A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  return svd.matrixV().col(3).hnormalized();
}

void ComputePoints3D(const std::vector<Track> &tracks,
                     const std::vector<Camera> &cams,
                     std::vector<Eigen::Vector3d> &tri_points) {

  tri_points.clear();
  for (std::vector<Track>::const_iterator citr = tracks.cbegin();
       citr != tracks.cend(); citr++) {
    if (1 < citr->size()) {
      tri_points.push_back(TriangulateTrack(*citr, cams));
    }
  }
}

Eigen::Matrix3d ComputeRotationMatrixViaEulerXYZ(double rx, double ry,
                                                 double rz) {
  return (Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ()))
      .toRotationMatrix();
}

std::vector<Eigen::Vector2d>
ProjectPoint(std::vector<Eigen::Vector3d> &points3d, Eigen::Matrix3d &R,
             Eigen::Vector3d &T, Eigen::Matrix3d &K) {

  cv::Mat cv_points3d(points3d.size(), 3, CV_32FC1);
  for (size_t idx = 0; idx < points3d.size(); idx++) {
    cv_points3d.at<float>(idx, 0) = points3d[idx](0);
    cv_points3d.at<float>(idx, 1) = points3d[idx](1);
    cv_points3d.at<float>(idx, 2) = points3d[idx](2);
    // LOG(INFO) << "x : " << points3d[idx](0) << ", y : " << points3d[idx](1)
    //          << ", z : " << points3d[idx](2);
  }

  cv::Mat cvR, cvT, cvK, rodrigues, cv_points2d;
  cv::eigen2cv(R, cvR);
  cv::eigen2cv(T, cvT);
  cv::eigen2cv(K, cvK);
  cv::Rodrigues(cvR, rodrigues);

  cv::projectPoints(cv_points3d, rodrigues, cvT, cvK, cv::noArray(),
                    cv_points2d);

  std::vector<Eigen::Vector2d> points2d;
  for (size_t idx = 0; idx < cv_points2d.rows; idx++) {
    Eigen::Vector2d p(cv_points2d.at<float>(idx, 0),
                      cv_points2d.at<float>(idx, 1));
    // LOG(INFO) << "x : " << p(0) << ", y : " << p(1);
    points2d.push_back(p);
  }
  return points2d;
}

} // namespace optimization