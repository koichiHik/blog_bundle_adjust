
// STL
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

// Google
#include <glog/logging.h>

// Eigen
#include <Eigen/Geometry>

// Original
#include "bundle_adjust.h"
#include "geometry.h"
#include "matrix_computation.h"
#include "rotations.h"

namespace {
std::pair<core::TrackId, core::CameraId> CreateKey(core::TrackId track_idx,
                                                   core::CameraId cam_idx) {
  return std::make_pair(track_idx, cam_idx);
}

template <typename T>
void CreatePointerVector(std::vector<T> &buff, std::vector<T *> &pointer_buff) {

  pointer_buff.clear();
  pointer_buff.resize(buff.size());
  for (size_t idx = 0; idx < buff.size(); idx++) {
    pointer_buff[idx] = &buff[idx];
  }
}

void DecomposeCameras(std::vector<core::Camera *> &cameras,
                      std::vector<Eigen::Matrix3d *> Ks,
                      std::vector<Eigen::Vector3d *> Ts,
                      std::vector<Eigen::Vector3d *> &Rots) {

  for (size_t cam_idx = 0; cam_idx < cameras.size(); cam_idx++) {

    Eigen::Matrix3d R;
    core::DecomposeCameraMatrix(*cameras[cam_idx], *Ks[cam_idx], R,
                                *Ts[cam_idx]);

    /*
    core::ComputeInternalCalibration(*cameras[cam_idx], *Ks[cam_idx], R);
    (*Ts[cam_idx]) = cameras[cam_idx]->block<3, 1>(0, 3);
    */

    // X. Convert to angle axis.
    Eigen::AngleAxisd axis(R);
    core::ConvertRotationMatrixToAngleAxis(R, *Rots[cam_idx]);
  }
}

double SumReprojectionError(
    const std::unordered_map<std::pair<core::TrackId, core::CameraId>,
                             Eigen::Vector2d *> &reproj_error) {

  double error_sum = 0;
  using Map = const std::unordered_map<std::pair<core::TrackId, core::CameraId>,
                                       Eigen::Vector2d *>;
  for (Map::const_iterator citr = reproj_error.cbegin();
       citr != reproj_error.cend(); citr++) {
    error_sum += citr->second->norm();
  }
  return error_sum;
}

} // namespace

namespace core {

using namespace std;

BundleAdjuster::BundleAdjuster(const std::vector<Camera> &cameras,
                               const std::vector<Track> &tracks,
                               const std::vector<Eigen::Vector3d> &points,
                               std::vector<Camera> &refined_cameras,
                               std::vector<Eigen::Vector3d> &refined_points) {

  // X. Copy input data.
  cam_num_ = cameras.size();
  track_num_ = tracks.size();
  buf_cameras_ = cameras;
  buf_tracks_ = tracks;
  buf_points_ = points;
  CreatePointerVector(buf_cameras_, cameras_);
  CreatePointerVector(buf_tracks_, tracks_);
  CreatePointerVector(buf_points_, points_);

  // X. Prepare cameras.
  buf_Ks_.resize(cam_num_, Eigen::Matrix3d::Zero());
  buf_Ts_.resize(cam_num_, Eigen::Vector3d::Zero());
  buf_Rots_.resize(cam_num_, Eigen::Vector3d::Zero());
  CreatePointerVector(buf_Ks_, Ks_);
  CreatePointerVector(buf_Ts_, Ts_);
  CreatePointerVector(buf_Rots_, Rots_);

  // X. Generate key set.
  CreateKeySet();

  // X. Initialize buffer.
  InitializeBuffer();

  // X. Run Levenberg Marquardt.
  RunOptimization(refined_cameras, refined_points);
}

void BundleAdjuster::CreateKeySet() {

  // X. Generate key set.
  for (size_t track_idx = 0; track_idx < tracks_.size(); track_idx++) {
    const Track &track = *tracks_[track_idx];
    for (Track::const_iterator citr = track.cbegin(); citr != track.cend();
         citr++) {
      size_t cam_idx = citr->first;
      key_set_.insert(make_pair(track_idx, cam_idx));

      // X. Track cam map.
      if (!track_cam_map_.count(track_idx)) {
        track_cam_map_.insert(
            make_pair(track_idx, std::unordered_set<CameraId>()));
      }
      track_cam_map_[track_idx].insert(cam_idx);

      // X. Cam track map.
      if (!cam_track_map_.count(cam_idx)) {
        cam_track_map_.insert(
            make_pair(cam_idx, std::unordered_set<TrackId>()));
      }
      cam_track_map_[cam_idx].insert(track_idx);
    }
  }
}

void BundleAdjuster::InitializeBuffer() {

  // X. Buffer initialize for (cam_idx, track_idx) matrix.
  {
    buf_reproj_errs_.resize(key_set_.size());
    buf_A_.resize(key_set_.size());
    buf_B_.resize(key_set_.size());
    buf_W_.resize(key_set_.size());
    buf_Y_.resize(key_set_.size());
    buf_Cov_.resize(key_set_.size(), MatrixCov::Identity());

    using Set = std::unordered_set<std::pair<TrackId, CameraId>>;
    size_t idx = 0;
    for (Set::const_iterator citr = key_set_.cbegin(); citr != key_set_.cend();
         citr++) {
      const std::pair<TrackId, CameraId> key = *citr;
      reproj_errs_.insert(make_pair(key, &buf_reproj_errs_[idx]));
      As_.insert(make_pair(key, &buf_A_[idx]));
      Bs_.insert(make_pair(key, &buf_B_[idx]));
      Ws_.insert(make_pair(key, &buf_W_[idx]));
      Ys_.insert(make_pair(key, &buf_Y_[idx]));
      CovInvs_.insert(make_pair(key, &buf_Cov_[idx]));
      idx++;
    }
  }

  // X. Buffer initalize for (cam_idx) matrix;
  {
    buf_U_.resize(cam_num_);
    buf_reproj_errs_on_A_.resize(cam_num_);
    buf_vecs_schur_.resize(cam_num_);
    size_t idx = 0;
    for (CameraId cam_idx = 0; cam_idx < cam_num_; cam_idx++) {
      Us_.push_back(&buf_U_[idx]);
      reproj_errs_on_A_.push_back(&buf_reproj_errs_on_A_[idx]);
      vecs_for_schur_.push_back(&buf_vecs_schur_[idx]);
      idx++;
    }
  }

  // X. Buffer initialize for (track_idx) matrix.
  {
    buf_V_.resize(track_num_);
    buf_reproj_errs_on_B_.resize(track_num_);
    size_t idx = 0;
    for (TrackId track_idx = 0; track_idx < track_num_; track_idx++) {
      Vs_.push_back(&buf_V_[idx]);
      reproj_errs_on_B_.push_back(&buf_reproj_errs_on_B_[idx]);
      idx++;
    }
  }

  // X. Buffer for linear system solver.
  {
    schur_comp_ = Eigen::MatrixXd(MatrixA::ColsAtCompileTime * cam_num_,
                                  MatrixA::ColsAtCompileTime * cam_num_);
    delta_ = Eigen::MatrixXd(MatrixA::ColsAtCompileTime * cam_num_ +
                                 MatrixB::ColsAtCompileTime * track_num_,
                             1);
    delta_A_ = Eigen::MatrixXd(MatrixA::ColsAtCompileTime * cam_num_, 1);
    delta_B_ = Eigen::MatrixXd(MatrixB::ColsAtCompileTime * track_num_, 1);
  }
}

void BundleAdjuster::RunOptimization(
    std::vector<Camera> &refined_cameras,
    std::vector<Eigen::Vector3d> &refined_points) {

  LOG(INFO) << "ComputeReprojectionError";
  ComputeReprojectionError(cameras_, tracks_, points_, reproj_errs_);

  LOG(INFO) << "Reprojection Error : " << SumReprojectionError(reproj_errs_);

  LOG(INFO) << "BuildLinearSystem";
  BuildLinearSystem();

  LOG(INFO) << "BuildLinearSystem";
  ProjectReprojectionErrorsOnParameters();

  LOG(INFO) << "SolveNormalEqutions";
  SolveNormalEqutions();

  LOG(INFO) << "UpdateParameters";
  UpdateParameters();

  LOG(INFO) << "ComputeReprojectionError";
  ComputeReprojectionError(cameras_, tracks_, points_, reproj_errs_);

  LOG(INFO) << "Reprojection Error : " << SumReprojectionError(reproj_errs_);

  LOG(INFO) << "Cycle Done";
}

void BundleAdjuster::ProjectReprojectionErrorsOnParameters() {

  // X. Compute reprojection on parameters.
  ComputeReprojectionErrorOnCameraParameters(As_, reproj_errs_, CovInvs_,
                                             reproj_errs_on_A_);

  ComputeReprojectionErrorOnPointParameters(Bs_, reproj_errs_, CovInvs_,
                                            reproj_errs_on_B_);
}

void BundleAdjuster::BuildLinearSystem() {

  DecomposeCameras(cameras_, Ks_, Ts_, Rots_);

  // X. Compute Matrix A, B, W
  {
    using Set = std::unordered_set<std::pair<TrackId, CameraId>>;
    for (Set::const_iterator citr = key_set_.begin(); citr != key_set_.end();
         citr++) {

      const std::pair<TrackId, CameraId> &key = *citr;
      MatrixA &A = *As_[key];
      MatrixB &B = *Bs_[key];
      MatrixW &W = *Ws_[key];
      ComputeMatrixA(*Ts_[key.second], *Rots_[key.second], *Ks_[key.second],
                     *points_[key.first], A);
      ComputeMatrixB(*Ts_[key.second], *Rots_[key.second], *Ks_[key.second],
                     *points_[key.first], B);
      ComputeMatrixW(A, B, *CovInvs_[key], W);
    }
  }

  // X. Compute U
  {
    using Map = std::unordered_map<CameraId, std::unordered_set<TrackId>>;
    for (Map::iterator itr = cam_track_map_.begin();
         itr != cam_track_map_.end(); itr++) {
      size_t cam_idx = itr->first;

      // X. Collect all A and Cov.
      std::vector<MatrixA *> tmp_As;
      std::vector<MatrixCov *> tmp_Covs;
      for (TrackId track_idx : itr->second) {
        Key key = make_pair(track_idx, cam_idx);
        tmp_As.push_back(As_[key]);
        tmp_Covs.push_back(CovInvs_[key]);
      }
      ComputeMatrixU(tmp_As, tmp_Covs, *Us_[cam_idx]);
    }
  }

  // X. Compute V
  {
    using Map = std::unordered_map<TrackId, std::unordered_set<CameraId>>;
    for (Map::const_iterator citr_t = track_cam_map_.cbegin();
         citr_t != track_cam_map_.cend(); citr_t++) {
      size_t track_idx = citr_t->first;

      // X. Collect
      std::vector<MatrixB *> tmp_Bs;
      std::vector<MatrixCov *> tmp_Covs;
      using Set = std::unordered_set<CameraId>;
      for (Set::const_iterator citr_c = citr_t->second.cbegin();
           citr_c != citr_t->second.cend(); citr_c++) {
        size_t cam_idx = *citr_c;
        Key key = make_pair(track_idx, cam_idx);
        tmp_Bs.push_back(Bs_[key]);
        tmp_Covs.push_back(CovInvs_[key]);
      }
      ComputeMatrixV(tmp_Bs, tmp_Covs, *Vs_[track_idx]);
    }
  }

  // X. Dumping factor.
  double myu = 1000000.0;
  AddDumpingFactor(myu, Us_, Vs_);
}

void BundleAdjuster::SolveNormalEqutions() {

  // X. Compute Y
  ComputeMatrixY(Ws_, Vs_, Ys_);

  // X. Compute Schur Complement.
  ComputeSchurComplement(Us_, Ys_, Ws_, schur_comp_);

  // X. Compute vector for Schur complement.
  ComputeVectorForSchur(Ys_, reproj_errs_on_A_, reproj_errs_on_B_,
                        vecs_for_schur_);

  // X. Solve for delta A.
  ComputeDeltaA(schur_comp_, vecs_for_schur_, delta_A_);

  std::cout << delta_A_;

  // X. Solve for delta B.
  ComputeDeltaB(Vs_, reproj_errs_on_B_, Ws_, delta_A_, delta_B_);

  // std::cout << delta_B_;

  // X. Merge delta A and delta B to delta.
  delta_.block(0, 0, MatrixA::ColsAtCompileTime * cam_num_, 1) = delta_A_;
  delta_.block(MatrixA::ColsAtCompileTime * cam_num_, 0,
               MatrixB::ColsAtCompileTime * track_num_, 1) = delta_B_;

  LOG(INFO) << "Delta norm : " << delta_.norm();
}

void BundleAdjuster::UpdateParameters() {

  for (size_t cam_idx = 0; cam_idx < cam_num_; cam_idx++) {
    Eigen::Block<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                 MatrixA::ColsAtCompileTime, 1>
        dA = delta_.block<MatrixA::ColsAtCompileTime, 1>(
            MatrixA::ColsAtCompileTime * cam_idx, 0);

    Eigen::Vector3d *T = Ts_[cam_idx];
    T->operator()(0) += dA(0, 0);
    T->operator()(1) += dA(1, 0);
    T->operator()(2) += dA(2, 0);

    Eigen::Vector3d *Rot = Rots_[cam_idx];
    Rot->operator()(0) += dA(3, 0);
    Rot->operator()(1) += dA(4, 0);
    Rot->operator()(2) += dA(5, 0);

    Eigen::Matrix3d *K = Ks_[cam_idx];
    K->operator()(0, 0) += dA(6, 0);
    K->operator()(1, 1) += dA(7, 0);
    K->operator()(0, 2) += dA(8, 0);
    K->operator()(1, 2) += dA(9, 0);
    K->operator()(0, 1) += dA(10, 0);

    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(*Rot, R);
    cameras_[cam_idx]->block<3, 3>(0, 0) = *K * R;
    cameras_[cam_idx]->block<3, 1>(0, 3) = *K * *T;
  }

  for (size_t track_idx = 0; track_idx < track_num_; track_idx++) {

    Eigen::Block<Eigen::Matrix<double, Eigen::Dynamic, 1>,
                 MatrixB::ColsAtCompileTime, 1>
        dB = delta_.block<MatrixB::ColsAtCompileTime, 1>(
            MatrixA::ColsAtCompileTime * cam_num_ +
                MatrixB::ColsAtCompileTime * track_idx,
            0);

    Eigen::Vector3d *p = points_[track_idx];
    p->operator()(0) += dB(0, 0);
    p->operator()(1) += dB(1, 0);
    p->operator()(2) += dB(2, 0);
  }
}

} // namespace core