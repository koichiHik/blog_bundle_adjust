

// STL
#include <iostream>

// Google
#include <glog/logging.h>

// Eigen
#include <Eigen/Dense>

// Original
#include "derivative_first_with_angle_axis.h"
#include "hash.h"
#include "matrix_computation.h"

namespace optimization {

// Camera Parameters, total 11
// tx, ty, tz, vx, vy, vz, fx, fy, cx, cy, s
//
// Matrix Aij : 2 x 11

void ComputeMatrixA(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    MatrixA &A) {

  Eigen::Vector3d x;
  {
    // X. Create rotation matrix.
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(rot, R);

    // X. Create camera matrix.
    Camera cam;
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T;
    cam = K * cam;
    x = cam * p.homogeneous();
  }

  // X. derivative for u.
  A(0, 0) = du_dtx(T, rot, K, p, x);
  A(0, 1) = du_dty(T, rot, K, p, x);
  A(0, 2) = du_dtz(T, rot, K, p, x);
  A(0, 3) = du_dvx(T, rot, K, p, x);
  A(0, 4) = du_dvy(T, rot, K, p, x);
  A(0, 5) = du_dvz(T, rot, K, p, x);
  A(0, 6) = du_dfx(T, rot, K, p, x);
  A(0, 7) = du_dfy(T, rot, K, p, x);
  A(0, 8) = du_dcx(T, rot, K, p, x);
  A(0, 9) = du_dcy(T, rot, K, p, x);
  A(0, 10) = du_ds(T, rot, K, p, x);

  // X. derivative for v.
  A(1, 0) = dv_dtx(T, rot, K, p, x);
  A(1, 1) = dv_dty(T, rot, K, p, x);
  A(1, 2) = dv_dtz(T, rot, K, p, x);
  A(1, 3) = dv_dvx(T, rot, K, p, x);
  A(1, 4) = dv_dvy(T, rot, K, p, x);
  A(1, 5) = dv_dvz(T, rot, K, p, x);
  A(1, 6) = dv_dfx(T, rot, K, p, x);
  A(1, 7) = dv_dfy(T, rot, K, p, x);
  A(1, 8) = dv_dcx(T, rot, K, p, x);
  A(1, 9) = dv_dcy(T, rot, K, p, x);
  A(1, 10) = dv_ds(T, rot, K, p, x);
}

// Point Parameters, total 3
// X, Y, Z
//
// Matrix Bij : 2 x 3

void ComputeMatrixB(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    MatrixB &B) {

  Eigen::Vector3d x;
  {
    // X. Create rotation matrix.
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(rot, R);

    // X. Create camera matrix.
    Camera cam;
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T;
    cam = K * cam;
    x = cam * p.homogeneous();
  }

  // X. derivative for u.
  B(0, 0) = du_dX(T, rot, K, p, x);
  B(0, 1) = du_dY(T, rot, K, p, x);
  B(0, 2) = du_dZ(T, rot, K, p, x);

  // X. derivative for v.
  B(1, 0) = dv_dX(T, rot, K, p, x);
  B(1, 1) = dv_dY(T, rot, K, p, x);
  B(1, 2) = dv_dZ(T, rot, K, p, x);
}

void ComputeMatrixW(const MatrixA &A, const MatrixB &B,
                    const Eigen::Matrix2d &Cov, MatrixW &W) {
  W = A.transpose() * Cov.inverse() * B;
}

void ComputeMatrixU(const std::vector<MatrixA *> &As,
                    const std::vector<Eigen::Matrix2d *> &Covs, MatrixU &U) {

  U.setZero();
  for (size_t idx = 0; idx < As.size(); idx++) {
    U = U + As[idx]->transpose() * Covs[idx]->inverse() * (*As[idx]);
  }
}

void ComputeMatrixV(const std::vector<MatrixB *> &Bs,
                    const std::vector<Eigen::Matrix2d *> &Cov, MatrixV &V) {

  V.setZero();
  for (size_t idx = 0; idx < Bs.size(); idx++) {
    V = V + Bs[idx]->transpose() * Cov[idx]->inverse() * (*Bs[idx]);
  }
}

void ComputeMatrixY(const std::unordered_map<Key, MatrixW *> &Ws,
                    const std::vector<MatrixV *> &Vs,
                    std::unordered_map<Key, MatrixY *> &Ys) {

  using Map = std::unordered_map<Key, MatrixY *>;
  for (Map::iterator itr = Ys.begin(); itr != Ys.end(); itr++) {
    size_t track_idx = itr->first.first;
    const MatrixV Vinv = Vs[track_idx]->inverse();
    (*Ys[itr->first]) = (*Ws.at(itr->first)) * Vinv;
  }
}

void ComputeVectorForSchur(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixY *> &Ys,
    const std::vector<VectorReprojOnA *> &reproj_error_on_As,
    const std::vector<VectorReprojOnB *> &reproj_error_on_Bs,
    std::vector<VectorForSchur *> &schurs_) {

  for (size_t cam_idx = 0; cam_idx < schurs_.size(); cam_idx++) {
    (*schurs_[cam_idx]) = (*reproj_error_on_As[cam_idx]);
  }

  using Map = std::unordered_map<std::pair<TrackId, CameraId>, MatrixY *>;
  for (Map::const_iterator citr = Ys.cbegin(); citr != Ys.cend(); citr++) {
    size_t track_idx = citr->first.first;
    size_t cam_idx = citr->first.second;
    (*schurs_[cam_idx]) -= (*citr->second) * (*reproj_error_on_Bs[track_idx]);
  }
}

void ComputeSchurComplement(
    const std::vector<MatrixU *> &Us,
    const std::unordered_map<Key, MatrixY *> &Ys,
    const std::unordered_map<Key, MatrixW *> &Ws,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &SchurComp) {

  SchurComp.setZero();

  size_t cam_num = Us.size();
  for (size_t row = 0; row < cam_num; row++) {
    for (size_t col = 0; col < cam_num; col++) {

      if (row == col) {
        SchurComp.block<MatrixA::ColsAtCompileTime, MatrixA::ColsAtCompileTime>(
            MatrixA::ColsAtCompileTime * row,
            MatrixA::ColsAtCompileTime * col) = *(Us[row]);
      }

      using Map = std::unordered_map<Key, MatrixY *>;
      for (Map::const_iterator citr = Ys.cbegin(); citr != Ys.cend(); citr++) {
        size_t cam_idx = citr->first.second;
        if (row == cam_idx) {
          SchurComp
              .block<MatrixA::ColsAtCompileTime, MatrixA::ColsAtCompileTime>(
                  MatrixA::ColsAtCompileTime * row,
                  MatrixA::ColsAtCompileTime * col) -=
              (*Ys.at(citr->first)) * (*Ws.at(citr->first)).transpose();
        }
      }
    }
  }
}

void AddDumpingFactor(double myu, std::vector<MatrixU *> &Us,
                      std::vector<MatrixV *> &Vs) {

  for (std::vector<MatrixU *>::iterator itr = Us.begin(); itr != Us.end();
       itr++) {
    for (int idx = 0; idx < MatrixU::RowsAtCompileTime; idx++) {
      (*itr)->operator()(idx, idx) += myu;
    }
  }

  for (std::vector<MatrixV *>::iterator itr = Vs.begin(); itr != Vs.end();
       itr++) {
    for (int idx = 0; idx < MatrixV::RowsAtCompileTime; idx++) {
      (*itr)->operator()(idx, idx) += myu;
    }
  }
}

void ComputeDeltaA(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &SchurComp,
    const std::vector<VectorForSchur *> &SchurVec,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaA) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> da(
      SchurVec.size() * VectorForSchur::RowsAtCompileTime, 1);

  for (size_t cam_idx = 0; cam_idx < SchurVec.size(); cam_idx++) {
    da.block<VectorForSchur::RowsAtCompileTime, 1>(
        VectorForSchur::RowsAtCompileTime * cam_idx, 0) = (*SchurVec[cam_idx]);
  }
  deltaA = SchurComp.colPivHouseholderQr().solve(da);
}

void ComputeDeltaB(const std::vector<MatrixV *> &Vs,
                   const std::vector<VectorReprojOnB *> &reproj_error_on_Bs,
                   const std::unordered_map<Key, MatrixW *> &Ws,
                   const Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaA,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaB) {

  deltaB.setZero();
  for (size_t track_idx = 0; track_idx < reproj_error_on_Bs.size();
       track_idx++) {
    deltaB.block<MatrixB::ColsAtCompileTime, 1>(
        MatrixB::ColsAtCompileTime * track_idx, 0) =
        *reproj_error_on_Bs[track_idx];
  }

  using Map = std::unordered_map<Key, MatrixW *>;
  for (Map::const_iterator citr = Ws.cbegin(); citr != Ws.cend(); citr++) {

    size_t track_idx = citr->first.first;
    size_t cam_idx = citr->first.second;

    deltaB.block<MatrixB::ColsAtCompileTime, 1>(
        MatrixB::ColsAtCompileTime * track_idx, 0) =
        (*Ws.at(citr->first)).transpose() *
        deltaA.block<MatrixA::ColsAtCompileTime, 1>(
            MatrixA::ColsAtCompileTime * cam_idx, 0);
  }

  for (size_t track_idx = 0; track_idx < Vs.size(); track_idx++) {
    deltaB.block<MatrixB::ColsAtCompileTime, 1>(
        MatrixB::ColsAtCompileTime * track_idx, 0) =
        (*Vs[track_idx]).inverse() *
        deltaB.block<MatrixB::ColsAtCompileTime, 1>(
            MatrixB::ColsAtCompileTime * track_idx, 0);
  }
}

/**
 *
 * Reprojection error computation.
 *
 */
void ComputeReprojectionError(
    const std::vector<Camera *> &cameras, const std::vector<Track *> &tracks,
    const std::vector<Eigen::Vector3d *> &points,
    std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error) {

  using Map =
      std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>;
  for (Map::iterator itr = reproj_error.begin(); itr != reproj_error.end();
       itr++) {
    size_t track_idx = itr->first.first;
    size_t cam_idx = itr->first.second;

    Eigen::Vector2d estimated =
        ((*cameras[cam_idx]) * (*points[track_idx]).homogeneous())
            .hnormalized();
    *(itr->second) = tracks[track_idx]->at(cam_idx) - estimated;
  }
}

void ComputeReprojectionErrorOnCameraParameters(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixA *> &As,
    const std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error,
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixCov *>
        &conv_invs,
    std::vector<VectorReprojOnA *> &reproj_errs_on_params) {

  using Vec = std::vector<VectorReprojOnA *>;
  for (Vec::iterator itr = reproj_errs_on_params.begin();
       itr != reproj_errs_on_params.end(); itr++) {
    (*itr)->setZero();
  }

  using Map = std::unordered_map<std::pair<TrackId, CameraId>, MatrixA *>;
  for (typename Map::const_iterator itr = As.cbegin(); itr != As.cend();
       itr++) {
    size_t cam_idx = itr->first.second;

    const MatrixA &A = *(itr->second);
    const Eigen::Vector2d &reproj_err = *(reproj_error.at(itr->first));
    const MatrixCov &cov_inv = *(conv_invs.at(itr->first));
    VectorReprojOnA &reproj_on_A = *(reproj_errs_on_params[cam_idx]);
    reproj_on_A = reproj_on_A + A.transpose() * cov_inv * reproj_err;
  }
}

void ComputeReprojectionErrorOnPointParameters(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixB *> &Bs,
    const std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error,
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixCov *>
        &conv_invs,
    std::vector<VectorReprojOnB *> &reproj_errs_on_params) {

  using Vec = std::vector<VectorReprojOnB *>;
  for (Vec::iterator itr = reproj_errs_on_params.begin();
       itr != reproj_errs_on_params.end(); itr++) {
    (*itr)->setZero();
  }

  using Map = std::unordered_map<std::pair<TrackId, CameraId>, MatrixB *>;
  for (typename Map::const_iterator itr = Bs.cbegin(); itr != Bs.cend();
       itr++) {
    size_t track_idx = itr->first.first;

    const MatrixB &B = *(itr->second);
    const Eigen::Vector2d &reproj_err = *(reproj_error.at(itr->first));
    const MatrixCov &cov_inv = *(conv_invs.at(itr->first));
    VectorReprojOnB &reproj_on_B = *(reproj_errs_on_params[track_idx]);
    reproj_on_B = reproj_on_B + B.transpose() * cov_inv * reproj_err;
  }
}

} // namespace optimization