

// STL
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Google
#include <glog/logging.h>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Original
#include "utility.h"

namespace utility {

void WritePoints3D(const std::string &path,
                   const std::vector<Eigen::Vector3d> &points3d) {

  std::ofstream write;
  write.open(path);
  for (std::vector<Eigen::Vector3d>::const_iterator citr = points3d.cbegin();
       citr != points3d.cend(); citr++) {
    write << std::fixed << std::setprecision(6) << (*citr)(0) << " "
          << (*citr)(1) << " " << (*citr)(2) << std::endl;
  }
}

void WritePoints2D(const std::string &path,
                   const std::vector<std::vector<Eigen::Vector2d>> &points2ds) {

  std::ofstream write;
  write.open(path);
  for (size_t track_idx = 0; track_idx < points2ds[0].size(); track_idx++) {
    for (size_t idx = 0; idx < points2ds.size(); idx++) {
      write << std::fixed << std::setprecision(6)
            << points2ds[idx][track_idx](0) << " "
            << points2ds[idx][track_idx](1);
      if (idx < points2ds.size() - 1) {
        write << " ";
      }
    }
    write << std::endl;
  }
}

void WriteCameraExtrinsicMatrix(const std::string &path,
                                const Eigen::Matrix3d &R,
                                const Eigen::Vector3d &T, bool overwrite) {

  std::ofstream write;
  if (overwrite) {
    write.open(path, std::ios::out);
  } else {
    write.open(path, std::ios::out | std::ios::app);
  }
  Eigen::Matrix<double, 3, 4> cam;
  cam.block<3, 3>(0, 0) = R;
  cam.block<3, 1>(0, 3) = T;
  write << std::fixed << std::setprecision(6) << cam(0, 0) << " " << cam(0, 1)
        << " " << cam(0, 2) << " " << cam(0, 3) << std::endl;
  write << std::fixed << std::setprecision(6) << cam(1, 0) << " " << cam(1, 1)
        << " " << cam(1, 2) << " " << cam(1, 3) << std::endl;
  write << std::fixed << std::setprecision(6) << cam(2, 0) << " " << cam(2, 1)
        << " " << cam(2, 2) << " " << cam(2, 3) << std::endl;
}

void WriteCameraIntrinsicMatrix(const std::string &path,
                                std::vector<size_t> cam_indices,
                                const Eigen::Matrix3d &K, bool overwrite) {

  std::ofstream write;
  if (overwrite) {
    write.open(path, std::ios::out);
  } else {
    write.open(path, std::ios::out | std::ios::app);
  }

  // X. camera intrinsic pair.
  for (size_t idx = 0; idx < cam_indices.size(); idx++) {
    write << cam_indices[idx];
    if (idx < cam_indices.size() - 1) {
      write << " ";
    }
  }
  write << std::endl;

  // X.
  write << std::fixed << std::setprecision(6) << K(0, 0) << " " << K(0, 1)
        << " " << K(0, 2) << std::endl;
  write << std::fixed << std::setprecision(6) << K(1, 0) << " " << K(1, 1)
        << " " << K(1, 2) << std::endl;
  write << std::fixed << std::setprecision(6) << K(2, 0) << " " << K(2, 1)
        << " " << K(2, 2) << std::endl;
}

std::string ReplaceString(const std::string &from, const std::string &to,
                          const std::string &input) {

  std::string replaced(input);
  std::string::size_type pos(replaced.find(from));

  while (pos != std::string::npos) {
    replaced.replace(pos, from.length(), to);
    pos = replaced.find(from, pos + to.length());
  }

  return replaced;
}

void ParseLineForDouble(const std::string &line, std::vector<double> &out) {
  std::string token;
  std::string tmp_line = ReplaceString("  ", " ", line);
  std::istringstream istr(tmp_line);

  while (std::getline(istr, token, ' ')) {

    try {
      if (token.size() > 0) {
        out.push_back(std::stod(token));
      }
    } catch (const std::invalid_argument &ia) {
      LOG(INFO) << "Invalid argument found, and skipped : " << token.size();
    }
  }
}

void ParseLineForSizeT(const std::string &line, std::vector<size_t> &out) {
  std::string token;
  std::string tmp_line = ReplaceString("  ", " ", line);
  std::istringstream istr(tmp_line);

  while (std::getline(istr, token, ' ')) {

    try {
      if (token.size() > 0) {
        out.push_back(std::stoul(token));
      }
    } catch (const std::invalid_argument &ia) {
      LOG(INFO) << "Invalid argument found, and skipped : " << token.size();
    }
  }
}

void LoadPoints3D(const std::string &path,
                  std::vector<Eigen::Vector3d> &points3d) {
  std::ifstream infile(path, std::ios::in);
  CHECK(infile.is_open()) << "File can not be opened.";

  std::string line;
  while (std::getline(infile, line)) {
    std::vector<double> row;
    ParseLineForDouble(line, row);
    Eigen::Vector3d p = Eigen::Map<Eigen::Vector3d>(row.data());
    points3d.push_back(p);
  }
}

void LoadCameraIntrinsicMatrix(const std::string &path,
                               std::map<size_t, size_t> &cam_intrinsic_map,
                               std::vector<Eigen::Matrix3d> &cam_intrinsics) {

  std::ifstream infile(path, std::ios::in);
  CHECK(infile.is_open()) << "File can not be opened.";

  std::string line;
  size_t intrinsic_idx = 0;
  while (std::getline(infile, line)) {

    // X. Load camera indices.
    std::vector<size_t> cam_indices;
    ParseLineForSizeT(line, cam_indices);
    for (size_t cam_idx : cam_indices) {
      cam_intrinsic_map[cam_idx] = intrinsic_idx;
    }
    intrinsic_idx++;

    // X. Load matrix.
    std::vector<double> row;
    std::getline(infile, line);
    ParseLineForDouble(line, row);
    std::getline(infile, line);
    ParseLineForDouble(line, row);
    std::getline(infile, line);
    ParseLineForDouble(line, row);

    Eigen::Matrix3d intrinsic =
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(row.data());

    double scale = 1.0;
    intrinsic(0, 0) /= scale;
    intrinsic(1, 1) /= scale;
    intrinsic(0, 2) /= scale;
    intrinsic(1, 2) /= scale;

    cam_intrinsics.push_back(intrinsic);
  }
}

void LoadCameraExtrinsicMatrix(
    const std::string &path,
    std::vector<Eigen::Matrix<double, 3, 4>> &cameras) {

  std::ifstream infile(path, std::ios::in);
  CHECK(infile.is_open()) << "File can not be opened.";

  std::string line;
  while (std::getline(infile, line)) {

    // X. 1st row.
    std::vector<double> row;
    ParseLineForDouble(line, row);
    std::getline(infile, line);
    ParseLineForDouble(line, row);
    std::getline(infile, line);
    ParseLineForDouble(line, row);

    Eigen::Matrix<double, 3, 4> cam =
        Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(row.data());

    cameras.push_back(cam);
  }
}

void LoadTracks(const std::string &path, int cam_no,
                std::vector<core::Track> &tracks) {

  std::ifstream infile(path, std::ios::in);
  CHECK(infile.is_open()) << "File can not be opened.";

  std::string line;
  while (std::getline(infile, line)) {

    std::vector<double> out;
    ParseLineForDouble(line, out);
    CHECK(out.size() == static_cast<size_t>(cam_no * 2))
        << "Read point number is not cam_no * 2 : " << cam_no * 2;

    double scale = 1.0;
    core::Track image_points;
    for (size_t cam_idx = 0; cam_idx < static_cast<size_t>(cam_no); cam_idx++) {
      if (out[cam_idx * 2] != -1) {
        image_points.insert(std::make_pair(
            cam_idx, Eigen::Vector2d(out[2 * cam_idx] / scale,
                                     out[2 * cam_idx + 1] / scale)));
      }
    }

    tracks.push_back(image_points);
  }
}

} // namespace utility