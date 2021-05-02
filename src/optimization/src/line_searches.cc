
// X. Google
#include <glog/logging.h>

// X. Original
#include <error_metric.h>
#include <optimizations.h>

//
#include "differential_operators.h"

namespace optimization {

double ComputeTrialValues(double alpha_low, double alpha_high) {

  double d_alpha = alpha_high - alpha_low;
  return alpha_low + 0.5 * d_alpha;
}

void Zoom(const std::vector<Track> &tracks_src,
          const std::vector<Eigen::Matrix3d> &K_src,
          const std::vector<Eigen::Vector3d> &T_src,
          const std::vector<Eigen::Vector3d> &Rot_src,
          const std::vector<Eigen::Vector3d> &points3d_src,
          const std::map<size_t, size_t> &extrinsic_intrinsic_map,
          double alpha_low, double alpha_high, double c1, double c2,
          size_t max_itr, double &alpha_star) {

  // X. Compute base value.
  double phi_0 = ComputeReprojectionError(
      tracks_src, K_src, extrinsic_intrinsic_map, Rot_src, T_src, points3d_src);
  Eigen::MatrixXd grad_0 = ComputeGradient(K_src, T_src, Rot_src, points3d_src,
                                           tracks_src, extrinsic_intrinsic_map);
  Eigen::MatrixXd step_dir = -grad_0.normalized();
  double der_phi_0 = grad_0.col(0).dot(step_dir.col(0));

  double phi_lo = ComputeReprojectionErrorWithStepLength(
      grad_0, -alpha_low, tracks_src, K_src, T_src, Rot_src, points3d_src,
      extrinsic_intrinsic_map);

  // X. Initial alpha
  double alpha_low_tmp_ = alpha_low;
  double alpha_high_tmp_ = alpha_high;

  size_t itr = 0;
  while (itr < max_itr) {
    // X. Find trial value.
    double alpha_j = ComputeTrialValues(alpha_low_tmp_, alpha_high_tmp_);

    // X. Update parameters based on alpha.
    std::vector<Eigen::Vector3d> T_tmp = T_src;
    std::vector<Eigen::Matrix3d> K_tmp = K_src;
    std::vector<Eigen::Vector3d> Rot_tmp = Rot_src;
    std::vector<Eigen::Vector3d> points3d_tmp = points3d_src;
    UpdateParameters(alpha_j * -grad_0, tracks_src, extrinsic_intrinsic_map,
                     K_tmp, T_tmp, Rot_tmp, points3d_tmp);
    double phi_j =
        ComputeReprojectionError(tracks_src, K_tmp, extrinsic_intrinsic_map,
                                 Rot_tmp, T_tmp, points3d_tmp);

    if (phi_j > phi_0 + c1 * alpha_j * der_phi_0 || phi_j >= phi_lo) {
      // X. Shrink upper window.
      alpha_high_tmp_ = alpha_j;
    } else {

      // X. Evaluate derivative value.
      Eigen::MatrixXd grad_j =
          ComputeGradient(K_tmp, T_tmp, Rot_tmp, points3d_tmp, tracks_src,
                          extrinsic_intrinsic_map);
      double der_phi_j = grad_j.col(0).dot(step_dir.col(0));

      // X. Find point satisfying wolfe condition.
      if (std::abs(der_phi_j) <= -c2 * der_phi_0) {
        alpha_star = alpha_j;
        break;
      }

      // X. Shrink upper window.
      if (der_phi_j * (alpha_high_tmp_ - alpha_low_tmp_) >= 0) {
        alpha_high_tmp_ = alpha_low_tmp_;
      }

      // X. Shrink lower window.
      alpha_low_tmp_ = alpha_j;
      phi_lo = phi_j;
    }
    itr++;
    if (max_itr <= itr) {
      LOG(FATAL) << "Zoom can not find target values.";
    }
  }
}

void LineSearchWithStrongWolfe(
    const std::vector<Track> &tracks_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map, size_t max_itr,
    double alpha_max, double c1, double c2, double &alpha_star) {

  // X. Temporary buffer.
  std::vector<Eigen::Vector3d> T_tmp = T_src;
  std::vector<Eigen::Matrix3d> K_tmp = K_src;
  std::vector<Eigen::Vector3d> Rot_tmp = Rot_src;
  std::vector<Eigen::Vector3d> points3d_tmp = points3d_src;

  // X. Compute base value.
  double phi_0 = ComputeReprojectionError(
      tracks_src, K_tmp, extrinsic_intrinsic_map, Rot_tmp, T_tmp, points3d_tmp);
  Eigen::MatrixXd grad_0 = ComputeGradient(K_tmp, T_tmp, Rot_tmp, points3d_tmp,
                                           tracks_src, extrinsic_intrinsic_map);
  Eigen::MatrixXd step_dir = -grad_0.normalized();
  double der_phi_0 = grad_0.col(0).dot(step_dir.col(0));

  // X. Start iteration.
  double phi_i = phi_0;
  double phi_i_1 = phi_0;
  double der_phi_i = der_phi_0;
  double der_phi_i_1 = der_phi_0;

  // X. Initial alpha.
  double alpha_i_1 = 0;
  double alpha_i = std::min(1.0, alpha_max);

  size_t itr = 0;
  while (itr < max_itr) {

    // X. Update parameters based on alpha.
    T_tmp = T_src;
    K_tmp = K_src;
    Rot_tmp = Rot_src;
    points3d_tmp = points3d_src;
    UpdateParameters(alpha_i * -grad_0, tracks_src, extrinsic_intrinsic_map,
                     K_tmp, T_tmp, Rot_tmp, points3d_tmp);

    // X. Compute reprojection error.
    phi_i = ComputeReprojectionError(tracks_src, K_tmp, extrinsic_intrinsic_map,
                                     Rot_tmp, T_tmp, points3d_tmp);

    // X. Condition 1. (Find range that target value exists.)
    if (phi_i > phi_0 + c1 * alpha_i * der_phi_0 ||
        (phi_i >= phi_i_1 && itr > 0)) {
      Zoom(tracks_src, K_src, T_src, Rot_src, points3d_src,
           extrinsic_intrinsic_map, alpha_i_1, alpha_i, c1, c2, max_itr,
           alpha_star);
      break;
    }

    // X. Compute gradient.
    Eigen::MatrixXd grad_i =
        ComputeGradient(K_tmp, T_tmp, Rot_tmp, points3d_tmp, tracks_src,
                        extrinsic_intrinsic_map);
    der_phi_i = grad_i.col(0).dot(step_dir.col(0));

    // X. Condition 2. (Find value that satisfied strong wolfe condition.)
    if (std::abs(der_phi_i) <= -c2 * der_phi_0) {
      alpha_star = alpha_i;
      break;
    }

    // X. Condition 3. (Find range that target value exists.)
    if (der_phi_i >= 0) {
      Zoom(tracks_src, K_src, T_src, Rot_src, points3d_src,
           extrinsic_intrinsic_map, alpha_i, alpha_i_1, c1, c2, max_itr,
           alpha_star);
      break;
    }

    // X. Bracket range can not be found. Retry with different alpha
    alpha_i_1 = alpha_i;
    alpha_i = std::min(8 * alpha_i, alpha_max);
    phi_i_1 = phi_i;
    der_phi_i_1 = der_phi_i;

    // X. Increment
    itr++;
    if (alpha_max < alpha_i || max_itr <= itr || alpha_i_1 == alpha_max) {
      LOG(FATAL) << "Line search algorithm did not converge.";
    }
  }
}

void LineSearchWithBackTracking(
    const std::vector<Track> &track_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map, size_t max_itr,
    double rho, double c, double &alpha_star) {

  CHECK(0 < rho && rho < 1) << "Value of rho has to be between 0 ~ 1";
  CHECK(0 < c && c <= 1) << "Value of c has to be between 0 ~ 1";

  // X. Compute base value.
  double phi_0 = ComputeReprojectionError(
      track_src, K_src, extrinsic_intrinsic_map, Rot_src, T_src, points3d_src);
  Eigen::MatrixXd grad_0 = ComputeGradient(K_src, T_src, Rot_src, points3d_src,
                                           track_src, extrinsic_intrinsic_map);
  Eigen::MatrixXd step_dir = -grad_0.normalized();
  double der_phi_0 = grad_0.col(0).dot(step_dir.col(0));
  // LOG(INFO) << "Gradient : " << der_phi_0;

  double alpha = 1.0;
  for (size_t itr = 0; itr < max_itr; itr++) {

    // X. Compute updated value.
    double phi_i = ComputeReprojectionErrorWithStepLength(
        grad_0, alpha, track_src, K_src, T_src, Rot_src, points3d_src,
        extrinsic_intrinsic_map);

    // LOG(INFO) << "Value : " << phi_i << ", " << phi_0 << ", alpha : " <<
    // alpha;

    if (phi_i <= phi_0 + c * alpha * der_phi_0) {
      break;
    }

    // X. Update contraction factor.
    alpha = rho * alpha;
  }
  alpha_star = alpha;
}

} // namespace optimization