
// STL
#include <vector>

// Google
#include <glog/logging.h>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

// Original
#include "common_def.h"
#include "derivative_first_with_angle_axis.h"
#include "derivative_u_second_with_angle_axis.h"
#include "derivative_v_second_with_angle_axis.h"
#include "differential_operators.h"

namespace optimization {

Eigen::MatrixXd MakeHessianPositiveDefiniteViaMultipleIdentity(
    const Eigen::MatrixXd& H, double multiple, double beta) {
  Eigen::VectorXd diagonals = H.diagonal();

  double tau = 0;
  double min_diag = diagonals.minCoeff();
  if (min_diag < 0) {
    tau = -min_diag + beta;
  }

  Eigen::LLT<Eigen::MatrixXd> cholesky;
  Eigen::MatrixXd H_tmp;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(H.rows(), H.cols());
  while (true) {
    H_tmp = H + tau * I;
    cholesky.compute(H_tmp);
    if (cholesky.info() == Eigen::ComputationInfo::Success) {
      break;
    }

    tau = std::max(multiple * tau, beta);
  }

  return H_tmp;
}

Eigen::MatrixXd ComputeHessian(
    const std::vector<Eigen::Matrix3d>& K,
    const std::vector<Eigen::Vector3d>& T,
    const std::vector<Eigen::Vector3d>& Rot,
    const std::vector<Eigen::Vector3d>& points3d,
    const std::vector<Track>& tracks,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map) {
  CHECK(points3d.size() == tracks.size())
      << "Number of points and tracks must be the same.";

  // X. Element count.
  size_t cam_ext_num = T.size();
  size_t cam_int_num = K.size();
  size_t track_num = tracks.size();
  size_t projected_pnt_num = 0;
  for (const Track& track : tracks) {
    projected_pnt_num += track.size();
  }

  // X. Prepare buffer.
  size_t param_num = cam_ext_num * kCameraExtrinsicParamsNum +
                     cam_int_num * kCameraIntrinsicParamsNum +
                     track_num * kPointParamsNum;
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(param_num, param_num);

  std::vector<Camera> cameras;
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {
    Camera cam;
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T[cam_ext_idx];
    cam = K[extrinsic_intrinsic_map.at(cam_ext_idx)] * cam;
    cameras.push_back(cam);
  }

  // X. Inner loop.
  for (size_t trk_idx = 0; trk_idx < track_num; trk_idx++) {
    for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {
      if (tracks[trk_idx].count(cam_ext_idx)) {
        size_t cam_int_idx = extrinsic_intrinsic_map.at(cam_ext_idx);

        Eigen::Vector2d meas = tracks[trk_idx].at(cam_ext_idx);
        Eigen::Vector3d x =
            cameras[cam_ext_idx] * points3d[trk_idx].homogeneous();
        Eigen::Vector2d proj = x.hnormalized();

        double u_res = (proj - meas)(0);
        double v_res = (proj - meas)(1);

        size_t row_idx_tx = cam_ext_idx * kCameraExtrinsicParamsNum;
        size_t row_idx_fx = cam_ext_num * kCameraExtrinsicParamsNum +
                            cam_int_idx * kCameraIntrinsicParamsNum;
        size_t row_idx_X = cam_ext_num * kCameraExtrinsicParamsNum +
                           cam_int_num * kCameraIntrinsicParamsNum +
                           trk_idx * kPointParamsNum;

        size_t col_idx_tx = cam_ext_idx * kCameraExtrinsicParamsNum;
        size_t col_idx_fx = cam_ext_num * kCameraExtrinsicParamsNum +
                            cam_int_idx * kCameraIntrinsicParamsNum;
        size_t col_idx_X = cam_ext_num * kCameraExtrinsicParamsNum +
                           cam_int_num * kCameraIntrinsicParamsNum +
                           trk_idx * kPointParamsNum;

        // 1st derivative x 1st derivative.
        {
          // tx.
          H(row_idx_tx, col_idx_tx) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 1) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 2) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 3) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 4) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx, col_idx_tx + 5) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx + 1) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx + 2) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx + 3) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx + 4) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx, col_idx_X) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx, col_idx_X + 1) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx, col_idx_X + 2) +=
              (du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // ty.
          H(row_idx_tx + 1, col_idx_tx) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_tx + 1) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_tx + 2) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_tx + 3) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 4) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 5) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_fx) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 1) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 2) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 3) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 4) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_X) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_X + 1) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_X + 2) +=
              (du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // tz.
          H(row_idx_tx + 2, col_idx_tx) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 1) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 2) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 3) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 4) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 5) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 2, col_idx_fx) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 2, col_idx_fx + 1) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 2) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 3) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 4) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx + 2, col_idx_X) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_X + 1) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_X + 2) +=
              (du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // vx.
          H(row_idx_tx + 3, col_idx_tx) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 1) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 2) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 3) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 4) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 5) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 3, col_idx_fx) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 1) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 2) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 3) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 4) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx + 3, col_idx_X) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_X + 1) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_X + 2) +=
              (du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // vy.
          H(row_idx_tx + 4, col_idx_tx) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 1) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 2) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 3) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 4) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 5) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 4, col_idx_fx) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 1) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 2) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 3) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 4) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx + 4, col_idx_X) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_X + 1) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_X + 2) +=
              (du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // vz.
          H(row_idx_tx + 5, col_idx_tx) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 1) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 2) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 3) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 4) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 5) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_tx + 5, col_idx_fx) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 1) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 2) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 3) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 4) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_tx + 5, col_idx_X) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_X + 1) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_X + 2) +=
              (du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // fx.
          H(row_idx_fx, col_idx_tx) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 1) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 2) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 3) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 4) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 5) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_fx, col_idx_fx) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 1) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 2) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 3) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 4) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_fx, col_idx_X) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_fx, col_idx_X + 1) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx, col_idx_X + 2) +=
              (du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // fy.
          H(row_idx_fx + 1, col_idx_tx) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 1) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 2) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 3) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 4) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 5) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_fx + 1, col_idx_fx) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 1) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 2) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 3) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 4) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_fx + 1, col_idx_X) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_X + 1) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_X + 2) +=
              (du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // cx.
          H(row_idx_fx + 2, col_idx_tx) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 1) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 2) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 3) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 4) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 5) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_fx + 2, col_idx_fx) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 1) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 2) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 3) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 4) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_fx + 2, col_idx_X) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_X + 1) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_X + 2) +=
              (du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // cy.
          H(row_idx_fx + 3, col_idx_tx) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 1) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 2) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 3) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 4) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 5) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_fx + 3, col_idx_fx) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 1) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 2) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 3) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 4) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_fx + 3, col_idx_X) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_X + 1) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_X + 2) +=
              (du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                      points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // s.
          H(row_idx_fx + 4, col_idx_tx) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 1) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 2) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 3) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 4) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 5) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_fx + 4, col_idx_fx) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 1) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 2) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 3) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 4) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X + 1) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X + 2) +=
              (du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // X.
          H(row_idx_X, col_idx_tx) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 1) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 2) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 3) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 4) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 5) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_X, col_idx_fx) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 1) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 2) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 3) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 4) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_X, col_idx_X) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_X, col_idx_X + 1) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_X, col_idx_X + 2) +=
              (du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // Y.
          H(row_idx_X + 1, col_idx_tx) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 1) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 2) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 3) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 4) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 5) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_X + 1, col_idx_fx) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 1) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 2) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 3) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 4) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_X + 1, col_idx_X) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_X + 1) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_X + 2) +=
              (du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          // Z.
          H(row_idx_X + 2, col_idx_tx) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 1) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 2) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 3) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 4) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 5) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));

          H(row_idx_X + 2, col_idx_fx) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 1) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 2) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 3) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                          points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 4) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));

          H(row_idx_X + 2, col_idx_X) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_X + 1) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_X + 2) +=
              (du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x) +
               dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x) *
                   dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                         points3d[trk_idx], x));
        }

        // 2nd derivatives.
        {
          // tx.
          H(row_idx_tx, col_idx_tx) +=
              u_res * (du2_dtxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 1) +=
              u_res * (du2_dtxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_tx + 2) +=
              u_res * (du2_dtxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx, col_idx_tx + 3) +=
              u_res * (du2_dtxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_tx + 4) +=
              u_res * (du2_dtxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_tx + 5) +=
              u_res * (du2_dtxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx, col_idx_fx) +=
              u_res * (du2_dtxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_fx + 1) +=
              u_res * (du2_dtxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_fx + 2) +=
              u_res * (du2_dtxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_fx + 3) +=
              u_res * (du2_dtxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_fx + 4) +=
              u_res * (du2_dtxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx, col_idx_X) +=
              u_res * (du2_dtxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_X + 1) +=
              u_res * (du2_dtxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx, col_idx_X + 2) +=
              u_res * (du2_dtxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // ty.
          H(row_idx_tx + 1, col_idx_tx) +=
              u_res * (du2_dtydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 1) +=
              u_res * (du2_dtydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 2) +=
              u_res * (du2_dtydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 3) +=
              u_res * (du2_dtydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 4) +=
              u_res * (du2_dtydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_tx + 5) +=
              u_res * (du2_dtydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_fx) +=
              u_res * (du2_dtydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 1) +=
              u_res * (du2_dtydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 2) +=
              u_res * (du2_dtydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 3) +=
              u_res * (du2_dtydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_fx + 4) +=
              u_res * (du2_dtyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 1, col_idx_X) +=
              u_res * (du2_dtydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_X + 1) +=
              u_res * (du2_dtydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 1, col_idx_X + 2) +=
              u_res * (du2_dtydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // tz.
          H(row_idx_tx + 2, col_idx_tx) +=
              u_res * (du2_dtzdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 1) +=
              u_res * (du2_dtzdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 2) +=
              u_res * (du2_dtzdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 3) +=
              u_res * (du2_dtzdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 4) +=
              u_res * (du2_dtzdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_tx + 5) +=
              u_res * (du2_dtzdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 2, col_idx_fx) +=
              u_res * (du2_dtzdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 1) +=
              u_res * (du2_dtzdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 2) +=
              u_res * (du2_dtzdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 3) +=
              u_res * (du2_dtzdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_fx + 4) +=
              u_res * (du2_dtzds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 2, col_idx_X) +=
              u_res * (du2_dtzdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_X + 1) +=
              u_res * (du2_dtzdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 2, col_idx_X + 2) +=
              u_res * (du2_dtzdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dtzdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // vx.
          H(row_idx_tx + 3, col_idx_tx) +=
              u_res * (du2_dvxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 1) +=
              u_res * (du2_dvxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 2) +=
              u_res * (du2_dvxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 3) +=
              u_res * (du2_dvxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 4) +=
              u_res * (du2_dvxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_tx + 5) +=
              u_res * (du2_dvxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 3, col_idx_fx) +=
              u_res * (du2_dvxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 1) +=
              u_res * (du2_dvxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 2) +=
              u_res * (du2_dvxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 3) +=
              u_res * (du2_dvxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_fx + 4) +=
              u_res * (du2_dvxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 3, col_idx_X) +=
              u_res * (du2_dvxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_X + 1) +=
              u_res * (du2_dvxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 3, col_idx_X + 2) +=
              u_res * (du2_dvxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // vy.
          H(row_idx_tx + 4, col_idx_tx) +=
              u_res * (du2_dvydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 1) +=
              u_res * (du2_dvydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 2) +=
              u_res * (du2_dvydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 3) +=
              u_res * (du2_dvydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 4) +=
              u_res * (du2_dvydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_tx + 5) +=
              u_res * (du2_dvydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 4, col_idx_fx) +=
              u_res * (du2_dvydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 1) +=
              u_res * (du2_dvydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 2) +=
              u_res * (du2_dvydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 3) +=
              u_res * (du2_dvydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_fx + 4) +=
              u_res * (du2_dvyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 4, col_idx_X) +=
              u_res * (du2_dvydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_X + 1) +=
              u_res * (du2_dvydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 4, col_idx_X + 2) +=
              u_res * (du2_dvydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // vz.
          H(row_idx_tx + 5, col_idx_tx) +=
              u_res * (du2_dvzdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 1) +=
              u_res * (du2_dvzdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 2) +=
              u_res * (du2_dvzdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 3) +=
              u_res * (du2_dvzdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 4) +=
              u_res * (du2_dvzdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_tx + 5) +=
              u_res * (du2_dvzdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 5, col_idx_fx) +=
              u_res * (du2_dvzdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 1) +=
              u_res * (du2_dvzdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 2) +=
              u_res * (du2_dvzdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 3) +=
              u_res * (du2_dvzdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_fx + 4) +=
              u_res * (du2_dvzds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_tx + 5, col_idx_X) +=
              u_res * (du2_dvzdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_X + 1) +=
              u_res * (du2_dvzdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_tx + 5, col_idx_X + 2) +=
              u_res * (du2_dvzdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dvzdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // fx.
          H(row_idx_fx, col_idx_tx) +=
              u_res * (du2_dfxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 1) +=
              u_res * (du2_dfxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 2) +=
              u_res * (du2_dfxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 3) +=
              u_res * (du2_dfxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 4) +=
              u_res * (du2_dfxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_tx + 5) +=
              u_res * (du2_dfxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx, col_idx_fx) +=
              u_res * (du2_dfxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 1) +=
              u_res * (du2_dfxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 2) +=
              u_res * (du2_dfxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 3) +=
              u_res * (du2_dfxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_fx + 4) +=
              u_res * (du2_dfxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx, col_idx_X) +=
              u_res * (du2_dfxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_X + 1) +=
              u_res * (du2_dfxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx, col_idx_X + 2) +=
              u_res * (du2_dfxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // fy.
          H(row_idx_fx + 1, col_idx_tx) +=
              u_res * (du2_dfydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 1) +=
              u_res * (du2_dfydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 2) +=
              u_res * (du2_dfydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 3) +=
              u_res * (du2_dfydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 4) +=
              u_res * (du2_dfydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_tx + 5) +=
              u_res * (du2_dfydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 1, col_idx_fx) +=
              u_res * (du2_dfydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 1) +=
              u_res * (du2_dfydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 2) +=
              u_res * (du2_dfydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 3) +=
              u_res * (du2_dfydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_fx + 4) +=
              u_res * (du2_dfyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 1, col_idx_X) +=
              u_res * (du2_dfydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_X + 1) +=
              u_res * (du2_dfydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 1, col_idx_X + 2) +=
              u_res * (du2_dfydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dfydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // cx.
          H(row_idx_fx + 2, col_idx_tx) +=
              u_res * (du2_dcxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 1) +=
              u_res * (du2_dcxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 2) +=
              u_res * (du2_dcxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 3) +=
              u_res * (du2_dcxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 4) +=
              u_res * (du2_dcxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_tx + 5) +=
              u_res * (du2_dcxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 2, col_idx_fx) +=
              u_res * (du2_dcxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 1) +=
              u_res * (du2_dcxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 2) +=
              u_res * (du2_dcxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 3) +=
              u_res * (du2_dcxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_fx + 4) +=
              u_res * (du2_dcxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 2, col_idx_X) +=
              u_res * (du2_dcxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_X + 1) +=
              u_res * (du2_dcxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 2, col_idx_X + 2) +=
              u_res * (du2_dcxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcxdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // cy.
          H(row_idx_fx + 3, col_idx_tx) +=
              u_res * (du2_dcydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 1) +=
              u_res * (du2_dcydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydty(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 2) +=
              u_res * (du2_dcydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 3) +=
              u_res * (du2_dcydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 4) +=
              u_res * (du2_dcydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_tx + 5) +=
              u_res * (du2_dcydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 3, col_idx_fx) +=
              u_res * (du2_dcydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 1) +=
              u_res * (du2_dcydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 2) +=
              u_res * (du2_dcydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 3) +=
              u_res * (du2_dcydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                  K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_fx + 4) +=
              u_res * (du2_dcyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcyds(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 3, col_idx_X) +=
              u_res * (du2_dcydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydX(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_X + 1) +=
              u_res * (du2_dcydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydY(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 3, col_idx_X + 2) +=
              u_res * (du2_dcydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dcydZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          // s.
          H(row_idx_fx + 4, col_idx_tx) +=
              u_res * (du2_dsdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 1) +=
              u_res * (du2_dsdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 2) +=
              u_res * (du2_dsdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 3) +=
              u_res * (du2_dsdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 4) +=
              u_res * (du2_dsdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_tx + 5) +=
              u_res * (du2_dsdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_fx + 4, col_idx_fx) +=
              u_res * (du2_dsdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 1) +=
              u_res * (du2_dsdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 2) +=
              u_res * (du2_dsdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 3) +=
              u_res * (du2_dsdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_fx + 4) +=
              u_res * (du2_dsds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X) +=
              u_res * (du2_dsdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X + 1) +=
              u_res * (du2_dsdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_fx + 4, col_idx_X + 2) +=
              u_res * (du2_dsdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dsdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          // X.
          H(row_idx_X, col_idx_tx) +=
              u_res * (du2_dXdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 1) +=
              u_res * (du2_dXdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 2) +=
              u_res * (du2_dXdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 3) +=
              u_res * (du2_dXdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 4) +=
              u_res * (du2_dXdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_tx + 5) +=
              u_res * (du2_dXdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X, col_idx_fx) +=
              u_res * (du2_dXdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 1) +=
              u_res * (du2_dXdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 2) +=
              u_res * (du2_dXdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 3) +=
              u_res * (du2_dXdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_fx + 4) +=
              u_res * (du2_dXds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X, col_idx_X) +=
              u_res * (du2_dXdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X, col_idx_X + 1) +=
              u_res * (du2_dXdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X, col_idx_X + 2) +=
              u_res * (du2_dXdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dXdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          // Y.
          H(row_idx_X + 1, col_idx_tx) +=
              u_res * (du2_dYdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 1) +=
              u_res * (du2_dYdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 2) +=
              u_res * (du2_dYdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 3) +=
              u_res * (du2_dYdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 4) +=
              u_res * (du2_dYdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_tx + 5) +=
              u_res * (du2_dYdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X + 1, col_idx_fx) +=
              u_res * (du2_dYdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 1) +=
              u_res * (du2_dYdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 2) +=
              u_res * (du2_dYdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 3) +=
              u_res * (du2_dYdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_fx + 4) +=
              u_res * (du2_dYds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X + 1, col_idx_X) +=
              u_res * (du2_dYdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_X + 1) +=
              u_res * (du2_dYdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 1, col_idx_X + 2) +=
              u_res * (du2_dYdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dYdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          // Z.
          H(row_idx_X + 2, col_idx_tx) +=
              u_res * (du2_dZdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 1) +=
              u_res * (du2_dZdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdty(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 2) +=
              u_res * (du2_dZdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdtz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 3) +=
              u_res * (du2_dZdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdvx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 4) +=
              u_res * (du2_dZdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdvy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_tx + 5) +=
              u_res * (du2_dZdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdvz(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X + 2, col_idx_fx) +=
              u_res * (du2_dZdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 1) +=
              u_res * (du2_dZdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 2) +=
              u_res * (du2_dZdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 3) +=
              u_res * (du2_dZdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                 K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_fx + 4) +=
              u_res * (du2_dZds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZds(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));

          H(row_idx_X + 2, col_idx_X) +=
              u_res * (du2_dZdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdX(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_X + 1) +=
              u_res * (du2_dZdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdY(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
          H(row_idx_X + 2, col_idx_X + 2) +=
              u_res * (du2_dZdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x)) +
              v_res * (dv2_dZdZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                K[cam_int_idx], points3d[trk_idx], x));
        }
      }
    }
  }

  return H;
}
}  // namespace optimization