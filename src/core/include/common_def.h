
#ifndef __COMMON_DEF_H__
#define __COMMON_DEF_H__

// STL
#include <map>

// Eigen
#include <Eigen/Core>

namespace optimization {

using TrackId = size_t;
using CameraId = size_t;
using Track = std::map<CameraId, Eigen::Vector2d>;
using Camera = Eigen::Matrix<double, 3, 4>;

static constexpr int kImageCoordNum = 2;
static constexpr int kCameraParamsNum = 11;
static constexpr int kPointParamsNum = 3;

using Key = std::pair<TrackId, CameraId>;
using MatrixCov = Eigen::Matrix<double, kImageCoordNum, kImageCoordNum>;
using MatrixA = Eigen::Matrix<double, kImageCoordNum, kCameraParamsNum>;
using MatrixB = Eigen::Matrix<double, kImageCoordNum, kPointParamsNum>;
using MatrixW = Eigen::Matrix<double, kCameraParamsNum, kPointParamsNum>;
using MatrixU = Eigen::Matrix<double, kCameraParamsNum, kCameraParamsNum>;
using MatrixV = Eigen::Matrix<double, kPointParamsNum, kPointParamsNum>;
using MatrixY = Eigen::Matrix<double, kCameraParamsNum, kPointParamsNum>;

using VectorReprojOnA = Eigen::Matrix<double, kCameraParamsNum, 1>;
using VectorReprojOnB = Eigen::Matrix<double, kPointParamsNum, 1>;
using VectorForSchur = VectorReprojOnA;

} // namespace optimization

#endif