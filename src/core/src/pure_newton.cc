
// STL
#include <iomanip>
#include <iostream>
#include <vector>

// Google
#include <glog/logging.h>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>

// Original
#include <error_metric.h>
#include <optimizations.h>

// Private header
#include "differential_operators.h"
#include "rotations.h"

namespace core {

void PureNewton::Optimize(
    const std::vector<Track> &tracks, const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Eigen::Vector3d> &points3d,
    std::vector<Camera> &extrinsics_dst,
    std::vector<Eigen::Matrix3d> &intrinsics_dst,
    std::vector<Eigen::Vector3d> &points3d_dst) {

  // X. Decompose camera matrix.
  std::vector<Eigen::Matrix3d> K = intrinsics;
  std::vector<Eigen::Vector3d> T(extrinsics.size(), Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> Rot(extrinsics.size(), Eigen::Vector3d::Zero());
  std::vector<Eigen::Vector3d> points3d_ = points3d;
  for (size_t cam_idx = 0; cam_idx < extrinsics.size(); cam_idx++) {
    T[cam_idx] = extrinsics[cam_idx].block<3, 1>(0, 3);
    core::ConvertRotationMatrixToAngleAxis(
        extrinsics[cam_idx].block<3, 3>(0, 0), Rot[cam_idx]);
  }

  // X. Start newton method.
  int loop_count = 0;
  double repj_new = ComputeReprojectionError(tracks, K, extrinsic_intrinsic_map,
                                             Rot, T, points3d_);
  double repj_old = repj_new;
  while (true) {

    // X. Compute gradient of reprojection error.
    Eigen::MatrixXd grad = core::ComputeGradient(K, T, Rot, points3d_, tracks,
                                                 extrinsic_intrinsic_map);
    Eigen::MatrixXd H = core::ComputeHessian(K, T, Rot, points3d_, tracks,
                                             extrinsic_intrinsic_map);

    // X. Make hessian positive definite.
    Eigen::MatrixXd Hp =
        core::MakeHessianPositiveDefiniteViaMultipleIdentity(H, 2, BETA);

    // X. Solve linear equations.
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(Hp);
    Eigen::VectorXd delta_params = lu_decomp.solve(-grad);
    UpdateParameters(delta_params, tracks, extrinsic_intrinsic_map, K, T, Rot,
                     points3d_);

    // X. Reprojection Error.
    repj_old = repj_new;
    repj_new = ComputeReprojectionError(tracks, K, extrinsic_intrinsic_map, Rot,
                                        T, points3d_);

    // X. Status report.
    if (true || delta_params.norm() < TERM_DELTA_NORM ||
        std::abs(repj_old - repj_new) <
            std::numeric_limits<double>::epsilon()) {

      double repj = ComputeReprojectionErrorWithStepLength(
          grad, 0.0, tracks, K, T, Rot, points3d_, extrinsic_intrinsic_map);
      LOG(INFO) << std::endl
                << "------------------ Current Status -------------------";
      LOG(INFO) << "Loop no            : " << loop_count;
      LOG(INFO) << "Delta param norm   : " << std::fixed << std::setprecision(9)
                << delta_params.norm() << " / " << TERM_DELTA_NORM;
      LOG(INFO) << "Reprojection Error : " << std::fixed
                << std::setprecision(20) << repj;

      if (delta_params.norm() < TERM_DELTA_NORM ||
          std::abs(repj_old - repj_new) < TERM_REPJ_DIFF) {
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
    core::ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);
    extrinsics_dst[cam_ext_idx].block<3, 3>(0, 0) = R;
    extrinsics_dst[cam_ext_idx].block<3, 1>(0, 3) = T[cam_ext_idx];
  }
}
} // namespace core