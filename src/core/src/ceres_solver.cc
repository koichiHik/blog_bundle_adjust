

// Google.
#include <ceres/ceres.h>
#include <glog/logging.h>

// X. Original
#include "ceres_reproj_error.h"
#include "rotations.h"
#include <optimizations.h>

void core::CeresSolver::Optimize(
    const std::vector<Track> &tracks, const std::vector<Camera> &extrinsics,
    const std::vector<Eigen::Matrix3d> &intrinsics,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map,
    const std::vector<Eigen::Vector3d> &points3d,
    std::vector<Camera> &extrinsics_dst,
    std::vector<Eigen::Matrix3d> &intrinsics_dst,
    std::vector<Eigen::Vector3d> &points3d_dst) {

  // X. Create buffer.
  std::vector<Eigen::Matrix<double, 5, 1>> K_tmp(intrinsics.size());
  std::vector<Eigen::Vector3d> points3d_tmp = points3d;
  std::vector<Eigen::Vector3d> Rot_tmp(extrinsics.size());
  std::vector<Eigen::Vector3d> T_tmp(extrinsics.size());

  // X. Convert intrinsic to vector
  for (size_t idx = 0; idx < intrinsics.size(); idx++) {
    const Eigen::Matrix3d &K = intrinsics[idx];
    K_tmp[idx](0) = K(0, 0);
    K_tmp[idx](1) = K(0, 1);
    K_tmp[idx](2) = K(0, 2);
    K_tmp[idx](3) = K(1, 1);
    K_tmp[idx](4) = K(1, 2);
  }

  // X. Convert camera to R and T.
  for (size_t idx = 0; idx < extrinsics.size(); idx++) {
    const Camera &cam = extrinsics[idx];
    core::ConvertRotationMatrixToAngleAxis(cam.block<3, 3>(0, 0), Rot_tmp[idx]);
    T_tmp[idx] = cam.block<3, 1>(0, 3);
  }

  // X. Create problem
  ceres::Problem problem;
  for (size_t cam_ext_idx = 0; cam_ext_idx < extrinsics.size(); cam_ext_idx++) {
    for (size_t trk_idx = 0; trk_idx < tracks.size(); trk_idx++) {

      if (tracks[trk_idx].count(cam_ext_idx)) {
        const Eigen::Vector2d meas = tracks[trk_idx].at(cam_ext_idx);
        ceres::CostFunction *cost_function =
            core::CeresReprojError::Create(meas(0), meas(1));

        Eigen::Vector3d &point = points3d_tmp[trk_idx];
        Eigen::Vector3d &Rot = Rot_tmp[cam_ext_idx];
        Eigen::Vector3d &T = T_tmp[cam_ext_idx];

        problem.AddResidualBlock(
            cost_function, nullptr, Rot_tmp[cam_ext_idx].data(),
            T_tmp[cam_ext_idx].data(),
            K_tmp[extrinsic_intrinsic_map.at(cam_ext_idx)].data(),
            points3d_tmp[trk_idx].data());
      }
    }
  }

  ceres::Solver::Options options;

  options.minimizer_type = ceres::MinimizerType::LINE_SEARCH;
  options.line_search_direction_type =
      ceres::LineSearchDirectionType::STEEPEST_DESCENT;

  options.linear_solver_type = ceres::DENSE_SCHUR;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << std::endl;

  // X. Convert extrinsic.
  extrinsics_dst.resize(extrinsics.size());
  for (size_t idx = 0; idx < T_tmp.size(); idx++) {
    Eigen::Matrix3d R_tmp;
    core::ConvertAngleAxisToRotationMatrix(Rot_tmp[idx], R_tmp);
    extrinsics_dst[idx].block<3, 3>(0, 0) = R_tmp;
    extrinsics_dst[idx].block<3, 1>(0, 3) = T_tmp[idx];
  }

  // X. Convert intrinsic to vector
  intrinsics_dst.resize(intrinsics.size());
  for (size_t idx = 0; idx < intrinsics.size(); idx++) {
    const Eigen::MatrixXd &K_vec = K_tmp[idx];
    Eigen::Matrix3d &K = intrinsics_dst[idx];

    K(0, 0) = K_vec(0);
    K(0, 1) = K_vec(1);
    K(0, 2) = K_vec(2);
    K(1, 0) = 0.0;
    K(1, 1) = K_vec(3);
    K(1, 2) = K_vec(4);
    K(2, 0) = 0.0;
    K(2, 1) = 0.0;
    K(2, 2) = 1.0;
  }

  // X. Points.
  points3d_dst = points3d_tmp;
}