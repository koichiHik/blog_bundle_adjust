

// STL
#include <iomanip>
#include <iostream>
#include <thread>

// Eigen
#include <Eigen/Dense>

// Google
#include <glog/logging.h>

// Original
#include <error_metric.h>
#include <geometry.h>
#include <optimizations.h>

// Private header
#include "differential_operators.h"
#include "rotations.h"

namespace {
void PrintStatus(const std::string& msg, int loop_count, double grad_norm,
                 double step_size, double repj_err, double term_grad_norm,
                 double term_step_length) {
  LOG(INFO) << std::endl << msg;
  LOG(INFO) << "Loop no            : " << loop_count;
  LOG(INFO) << "Gradient norm      : " << grad_norm << " / " << term_grad_norm;
  LOG(INFO) << "Step               : " << step_size << " / "
            << term_step_length;
  LOG(INFO) << "Reprojection Error : " << repj_err;
}
}  // namespace

namespace optimization {

void GradientDescent::Optimize(
    const std::vector<Track>& tracks_src,
    const std::vector<Camera>& extrinsics_src,
    const std::vector<Eigen::Matrix3d>& intrinsics_src,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map_src,
    const std::vector<Eigen::Vector3d>& points3d_src,
    std::vector<Camera>& extrinsics_dst,
    std::vector<Eigen::Matrix3d>& intrinsics_dst,
    std::vector<Eigen::Vector3d>& points3d_dst) {
  // X. Decompose camera matrix.
  std::vector<Eigen::Matrix3d> K = intrinsics_src;
  std::vector<Eigen::Vector3d> T(extrinsics_src.size(),
                                 Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> Rot(extrinsics_src.size(),
                                   Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> points3d_ = points3d_src;
  for (size_t cam_idx = 0; cam_idx < extrinsics_src.size(); cam_idx++) {
    T[cam_idx] = extrinsics_src[cam_idx].block<3, 1>(0, 3);
    optimization::ConvertRotationMatrixToAngleAxis(
        extrinsics_src[cam_idx].block<3, 3>(0, 0), Rot[cam_idx]);
  }

  // X. Start gradient descent.
  int loop_count = 0;
  double step = INITIAL_STEP;
  double last_repj_err;
  while (true) {
    // X. Compute gradient of reprojection error.
    Eigen::MatrixXd grad = optimization::ComputeGradient(
        K, T, Rot, points3d_, tracks_src, extrinsic_intrinsic_map_src);

    // X. Compute New and Old reprojection error.
    double repj_old = optimization::ComputeReprojectionErrorWithStepLength(
        -grad, 0, tracks_src, K, T, Rot, points3d_,
        extrinsic_intrinsic_map_src);
    double repj_new = optimization::ComputeReprojectionErrorWithStepLength(
        -grad, step, tracks_src, K, T, Rot, points3d_,
        extrinsic_intrinsic_map_src);

    if (loop_count % 1000 == 0) {
      last_repj_err = repj_old;
      std::string msg = "------------------ Current Status -------------------";
      PrintStatus(msg, loop_count, grad.norm(), step, repj_old,
                  TERM_GRADIENT_NORM, TERM_STEP_LENGTH);
    }

    if (repj_new < repj_old) {
      while (true) {
        // X. Update old parameters.
        UpdateParameters(step * (-grad), tracks_src,
                         extrinsic_intrinsic_map_src, K, T, Rot, points3d_);
        repj_old = repj_new;

        // X. Compute new reprojection error.
        step *= 2.0;
        repj_new = optimization::ComputeReprojectionErrorWithStepLength(
            -grad, step, tracks_src, K, T, Rot, points3d_,
            extrinsic_intrinsic_map_src);

        // X. Break if passes minimum.
        if (repj_old <= repj_new) {
          step *= 0.5;
          break;
        }
      }
    } else {
      while (true) {
        // X. Compute new reprojection error.
        step *= 0.5;
        repj_new = optimization::ComputeReprojectionErrorWithStepLength(
            -grad, step, tracks_src, K, T, Rot, points3d_,
            extrinsic_intrinsic_map_src);

        // X. Break if finds smaller value.
        if (repj_new <= repj_old) {
          UpdateParameters(step * (-grad), tracks_src,
                           extrinsic_intrinsic_map_src, K, T, Rot, points3d_);
          step *= 2.0;
          break;
        }
      }
    }

    // X. Check termination condition.
    if (grad.norm() < TERM_GRADIENT_NORM || step < TERM_STEP_LENGTH ||
        std::abs(last_repj_err - repj_new) < TERM_REPJ_DIFF) {
      // X. Final status report.
      std::string msg = "------------------ Final Status -------------------";
      PrintStatus(msg, loop_count, grad.norm(), step, repj_old,
                  TERM_GRADIENT_NORM, TERM_STEP_LENGTH);

      break;
    }

    loop_count++;
  }

  intrinsics_dst = K;
  points3d_dst = points3d_;
  extrinsics_dst.resize(T.size());
  for (size_t cam_ext_idx = 0; cam_ext_idx < T.size(); cam_ext_idx++) {
    Eigen::Matrix3d R;
    optimization::ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);
    extrinsics_dst[cam_ext_idx].block<3, 3>(0, 0) = R;
    extrinsics_dst[cam_ext_idx].block<3, 1>(0, 3) = T[cam_ext_idx];
  }
}

void GradientDescentWithLineSearch::Optimize(
    const std::vector<Track>& tracks_src,
    const std::vector<Camera>& extrinsics_src,
    const std::vector<Eigen::Matrix3d>& intrinsics_src,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map_src,
    const std::vector<Eigen::Vector3d>& points3d_src,
    std::vector<Camera>& extrinsics_dst,
    std::vector<Eigen::Matrix3d>& intrinsics_dst,
    std::vector<Eigen::Vector3d>& points3d_dst) {
  // X. Decompose camera matrix.
  std::vector<Eigen::Matrix3d> K = intrinsics_src;
  std::vector<Eigen::Vector3d> T(extrinsics_src.size(),
                                 Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> Rot(extrinsics_src.size(),
                                   Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> points3d_ = points3d_src;
  for (size_t cam_idx = 0; cam_idx < extrinsics_src.size(); cam_idx++) {
    T[cam_idx] = extrinsics_src[cam_idx].block<3, 1>(0, 3);
    optimization::ConvertRotationMatrixToAngleAxis(
        extrinsics_src[cam_idx].block<3, 3>(0, 0), Rot[cam_idx]);
  }

  // X. Start gradient descent.
  int loop_count = 0;
  double alpha = 0;
  while (true) {
    // X. Solve line search for the gradient direction.
    bool line_search_success = LineSearchWithStrongWolfe(
        tracks_src, K, T, Rot, points3d_, extrinsic_intrinsic_map_src,
        LINE_SEARCH_MAX_ITR, LINE_SEARCH_MIN_STEP, LINE_SEARCH_ALPHA_MAX,
        LINE_SEARCH_C1, LINE_SEARCH_C2, alpha);
    if (!line_search_success) {
      LOG(INFO) << "Line search did not converge. Optimization failed.";
      break;
    }

    // X. Reflect result from line search.
    Eigen::MatrixXd grad = optimization::ComputeGradient(
        K, T, Rot, points3d_, tracks_src, extrinsic_intrinsic_map_src);
    UpdateParameters(alpha * -grad, tracks_src, extrinsic_intrinsic_map_src, K,
                     T, Rot, points3d_);

    // X. Compute new gradient.
    Eigen::MatrixXd new_grad = optimization::ComputeGradient(
        K, T, Rot, points3d_, tracks_src, extrinsic_intrinsic_map_src);

    if (loop_count % 1000 == 0 || new_grad.norm() < TERM_GRADIENT_NORM) {
      // X. Compute New and Old reprojection error.
      double repj = optimization::ComputeReprojectionError(
          tracks_src, K, extrinsic_intrinsic_map_src, Rot, T, points3d_);

      if (TERM_GRADIENT_NORM < new_grad.norm()) {
        std::string msg =
            "------------------ Current Status -------------------";
        PrintStatus(msg, loop_count, new_grad.norm(), alpha, repj,
                    TERM_GRADIENT_NORM, TERM_STEP_LENGTH);
      } else {
        // X. Final status report.
        std::string msg = "------------------ Final Status -------------------";
        PrintStatus(msg, loop_count, new_grad.norm(), alpha, repj,
                    TERM_GRADIENT_NORM, TERM_STEP_LENGTH);
        break;
      }
    }

    loop_count++;
  }

  intrinsics_dst = K;
  points3d_dst = points3d_;
  extrinsics_dst.resize(T.size());
  for (size_t cam_ext_idx = 0; cam_ext_idx < T.size(); cam_ext_idx++) {
    Eigen::Matrix3d R;
    optimization::ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);
    extrinsics_dst[cam_ext_idx].block<3, 3>(0, 0) = R;
    extrinsics_dst[cam_ext_idx].block<3, 1>(0, 3) = T[cam_ext_idx];
  }
}

}  // namespace optimization