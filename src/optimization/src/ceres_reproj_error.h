
#ifndef __CERES_REPROJ_ERROR_H__
#define __CERES_REPROJ_ERROR_H__

// Ceres
#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/rotation.h>

// Original

namespace optimization {

struct CeresReprojError {

  CeresReprojError(double observed_x, double observed_y)
      : observed_x_(observed_x), observed_y_(observed_y) {}

  template <typename T>
  bool operator()(const T *const rot, const T *const trans, const T *const K,
                  const T *const point, T *residuals) const {

    // X. Convert point to camera coordinate system.
    T p[3];
    ceres::AngleAxisRotatePoint(rot, point, p);
    p[0] += trans[0];
    p[1] += trans[1];
    p[2] += trans[2];

    // X. Predicted image coord.
    // K[0] = fx, K[1] = s, K[2] = cx, K[3] = fy, K[4] = cy
    T xp = (K[0] * p[0] + K[1] * p[1] + K[2] * p[2]) / p[2];
    T yp = (K[3] * p[1] + K[4] * p[2]) / p[2];

    residuals[0] = xp - T(observed_x_);
    residuals[1] = yp - T(observed_y_);

    return true;
  }

  static ceres::CostFunction *Create(const double observed_x,
                                     const double observed_y) {
    return (new ceres::AutoDiffCostFunction<CeresReprojError, 2, 3, 3, 5, 3>(
        new CeresReprojError(observed_x, observed_y)));
  }

  double observed_x_;
  double observed_y_;
};

} // namespace optimization

#endif