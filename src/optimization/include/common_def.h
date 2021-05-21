
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

static constexpr size_t kCameraExtrinsicParamsNum = 6;
static constexpr size_t kCameraIntrinsicParamsNum = 5;
static constexpr size_t kPointParamsNum = 3;

}  // namespace optimization

#endif