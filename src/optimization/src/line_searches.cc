
// X. Google
#include <glog/logging.h>

// X. Original
#include <error_metric.h>
#include <optimizations.h>

//
#include "differential_operators.h"

namespace optimization {

double ComputeTrialValues(double alpha_low, double alpha_high) {
  return (alpha_low + alpha_high) / 2.0;
}

bool Zoom(const std::vector<Track>& tracks_src,
          const std::vector<Eigen::Matrix3d>& K_src,
          const std::vector<Eigen::Vector3d>& T_src,
          const std::vector<Eigen::Vector3d>& Rot_src,
          const std::vector<Eigen::Vector3d>& points3d_src,
          const std::map<size_t, size_t>& extrinsic_intrinsic_map,
          const Eigen::MatrixXd& gauge_fix_matrix, double alpha_low,
          double alpha_high, double c1, double c2, size_t max_itr,
          double min_step, double& alpha_star) {
  // X. Compute base value & base gradient.
  double phi_0 = ComputeReprojectionError(
      tracks_src, K_src, extrinsic_intrinsic_map, Rot_src, T_src, points3d_src);
  Eigen::MatrixXd grad_0 = gauge_fix_matrix.transpose() *
                           ComputeGradient(K_src, T_src, Rot_src, points3d_src,
                                           tracks_src, extrinsic_intrinsic_map);
  double grad_0_max_norm = grad_0.col(0).lpNorm<Eigen::Infinity>();

  // X. Compute directional derivative.
  Eigen::MatrixXd step_dir = -grad_0.normalized();
  double der_phi_0 = grad_0.col(0).dot(step_dir.col(0));

  // X. Compute value at alpha_low.
  double phi_lo = ComputeReprojectionErrorWithStepLength(
      -gauge_fix_matrix * grad_0, alpha_low, tracks_src, K_src, T_src, Rot_src,
      points3d_src, extrinsic_intrinsic_map);

  // X. Initial alpha
  double alpha_low_tmp_ = alpha_low;
  double alpha_high_tmp_ = alpha_high;

  // X. Result
  bool success = false;

  for (size_t itr = 0; itr < max_itr; itr++) {
    // X. Bracket Size Check
    if ((std::abs(alpha_high_tmp_ - alpha_low_tmp_) * grad_0_max_norm <
         min_step)) {
      success = false;
      alpha_star = alpha_low;
      break;
    }

    // X. Find trial value.
    double alpha_j = ComputeTrialValues(alpha_low_tmp_, alpha_high_tmp_);

    // X. Update parameters based on alpha.
    std::vector<Eigen::Vector3d> T_tmp = T_src;
    std::vector<Eigen::Matrix3d> K_tmp = K_src;
    std::vector<Eigen::Vector3d> Rot_tmp = Rot_src;
    std::vector<Eigen::Vector3d> points3d_tmp = points3d_src;
    UpdateParameters(alpha_j * gauge_fix_matrix * -grad_0, tracks_src,
                     extrinsic_intrinsic_map, K_tmp, T_tmp, Rot_tmp,
                     points3d_tmp);
    double phi_j =
        ComputeReprojectionError(tracks_src, K_tmp, extrinsic_intrinsic_map,
                                 Rot_tmp, T_tmp, points3d_tmp);

    if (phi_j > phi_0 + c1 * alpha_j * der_phi_0 || phi_j >= phi_lo) {
      // X. Shrink upper window.
      alpha_high_tmp_ = alpha_j;
    } else {
      // X. Evaluate derivative value.
      Eigen::MatrixXd grad_j =
          gauge_fix_matrix.transpose() *
          ComputeGradient(K_tmp, T_tmp, Rot_tmp, points3d_tmp, tracks_src,
                          extrinsic_intrinsic_map);
      double der_phi_j = grad_j.col(0).dot(step_dir.col(0));

      // X. Find point satisfying wolfe condition.
      if (std::abs(der_phi_j) <= -c2 * der_phi_0) {
        alpha_star = alpha_j;
        success = true;
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
  }

  return success;
}

bool LineSearchWithStrongWolfe(
    const std::vector<Track>& tracks_src,
    const std::vector<Eigen::Matrix3d>& K_src,
    const std::vector<Eigen::Vector3d>& T_src,
    const std::vector<Eigen::Vector3d>& Rot_src,
    const std::vector<Eigen::Vector3d>& points3d_src,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map,
    const Eigen::MatrixXd& gauge_fix_matrix, size_t max_itr, double min_step,
    double alpha_max, double c1, double c2, double& alpha_star) {
  // X. Temporary buffer.
  std::vector<Eigen::Vector3d> T_tmp = T_src;
  std::vector<Eigen::Matrix3d> K_tmp = K_src;
  std::vector<Eigen::Vector3d> Rot_tmp = Rot_src;
  std::vector<Eigen::Vector3d> points3d_tmp = points3d_src;

  // X. Compute base value & base gradient.
  double phi_0 = ComputeReprojectionError(
      tracks_src, K_tmp, extrinsic_intrinsic_map, Rot_tmp, T_tmp, points3d_tmp);
  Eigen::MatrixXd grad_0 = gauge_fix_matrix.transpose() *
                           ComputeGradient(K_tmp, T_tmp, Rot_tmp, points3d_tmp,
                                           tracks_src, extrinsic_intrinsic_map);
  double grad_0_max_norm = grad_0.col(0).lpNorm<Eigen::Infinity>();

  // X. Compute directional derivative.
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

  // X. Status
  bool success = false;

  for (size_t itr = 0; itr < max_itr; itr++) {
    // X. Copy from original..
    T_tmp = T_src;
    K_tmp = K_src;
    Rot_tmp = Rot_src;
    points3d_tmp = points3d_src;

    // X. Update parameters based on alpha.
    UpdateParameters(alpha_i * gauge_fix_matrix * -grad_0, tracks_src,
                     extrinsic_intrinsic_map, K_tmp, T_tmp, Rot_tmp,
                     points3d_tmp);

    // X. Compute function value at alpha_i.
    phi_i = ComputeReprojectionError(tracks_src, K_tmp, extrinsic_intrinsic_map,
                                     Rot_tmp, T_tmp, points3d_tmp);

    // X. Condition 1.
    if (phi_i > phi_0 + c1 * alpha_i * der_phi_0 ||
        (phi_i >= phi_i_1 && itr > 0)) {
      success = Zoom(tracks_src, K_src, T_src, Rot_src, points3d_src,
                     extrinsic_intrinsic_map, gauge_fix_matrix, alpha_i_1,
                     alpha_i, c1, c2, max_itr, min_step, alpha_star);
      break;
    }

    // X. Compute gradient at alpha_i.
    Eigen::MatrixXd grad_i =
        gauge_fix_matrix.transpose() * ComputeGradient(K_tmp, T_tmp, Rot_tmp,
                                                       points3d_tmp, tracks_src,
                                                       extrinsic_intrinsic_map);
    der_phi_i = grad_i.col(0).dot(step_dir.col(0));

    // X. Condition 2.
    if (std::abs(der_phi_i) <= -c2 * der_phi_0) {
      alpha_star = alpha_i;
      success = true;
      break;
    }

    // X. Condition 3.
    if (der_phi_i >= 0) {
      success = Zoom(tracks_src, K_src, T_src, Rot_src, points3d_src,
                     extrinsic_intrinsic_map, gauge_fix_matrix, alpha_i,
                     alpha_i_1, c1, c2, max_itr, min_step, alpha_star);
      break;
    }

    // X. Bracket range can not be found. Retry with different alpha
    alpha_i_1 = alpha_i;
    alpha_i = std::min(2 * alpha_i, alpha_max);
    phi_i_1 = phi_i;
    der_phi_i_1 = der_phi_i;

    // X. Bracket Size Check
    if ((std::abs(alpha_i - alpha_i_1) * grad_0_max_norm < min_step) ||
        (alpha_i_1 == alpha_max)) {
      success = false;
      alpha_star = alpha_i_1;
      break;
    }
  }

  return success;
}

void LineSearchWithBackTracking(
    const std::vector<Track>& track_src,
    const std::vector<Eigen::Matrix3d>& K_src,
    const std::vector<Eigen::Vector3d>& T_src,
    const std::vector<Eigen::Vector3d>& Rot_src,
    const std::vector<Eigen::Vector3d>& points3d_src,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map, size_t max_itr,
    double rho, double c, double& alpha_star) {
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

}  // namespace optimization