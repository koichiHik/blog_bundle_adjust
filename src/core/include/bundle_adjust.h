
#ifndef __BUNDLE_ADJUST_H__
#define __BUNDLE_ADJUST_H__

// STL
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Eigen
#include <Eigen/Core>

// Original
#include "hash.h"
#include <common_def.h>

namespace core {

class BundleAdjuster {
public:
  BundleAdjuster(const std::vector<Camera> &cameras,
                 const std::vector<Track> &tracks,
                 const std::vector<Eigen::Vector3d> &points,
                 std::vector<Camera> &refined_cameras,
                 std::vector<Eigen::Vector3d> &refined_points);

private:
  void CreateKeySet();

  void InitializeBuffer();

  void RunOptimization(std::vector<Camera> &refined_cameras,
                       std::vector<Eigen::Vector3d> &refined_points);

  void ProjectReprojectionErrorsOnParameters();

  void BuildLinearSystem();

  void SolveNormalEqutions();

  void UpdateParameters();

private:
  //
  size_t cam_num_;
  size_t track_num_;

  // X. Key related.
  std::unordered_set<std::pair<TrackId, CameraId>> key_set_;
  std::unordered_map<CameraId, std::unordered_set<TrackId>> cam_track_map_;
  std::unordered_map<TrackId, std::unordered_set<CameraId>> track_cam_map_;

  // X. Buffer for input data.
  std::vector<Camera *> cameras_;
  std::vector<Camera> buf_cameras_;
  std::vector<Track *> tracks_;
  std::vector<Track> buf_tracks_;
  std::vector<Eigen::Vector3d *> points_;
  std::vector<Eigen::Vector3d> buf_points_;
  std::unordered_map<Key, MatrixCov *> CovInvs_;
  std::vector<MatrixCov> buf_Cov_;

  // X. Buffer for cameras.
  std::vector<Eigen::Matrix3d *> Ks_;
  std::vector<Eigen::Matrix3d> buf_Ks_;
  std::vector<Eigen::Vector3d *> Ts_;
  std::vector<Eigen::Vector3d> buf_Ts_;
  std::vector<Eigen::Vector3d *> Rots_;
  std::vector<Eigen::Vector3d> buf_Rots_;

  // X. Matrix Reference.
  std::unordered_map<Key, MatrixA *> As_;
  std::vector<MatrixA> buf_A_;
  std::unordered_map<Key, MatrixB *> Bs_;
  std::vector<MatrixB> buf_B_;
  std::unordered_map<Key, MatrixW *> Ws_;
  std::vector<MatrixW> buf_W_;
  std::unordered_map<Key, MatrixY *> Ys_;
  std::vector<MatrixY> buf_Y_;
  std::vector<MatrixU *> Us_;
  std::vector<MatrixU> buf_U_;
  std::vector<MatrixV *> Vs_;
  std::vector<MatrixV> buf_V_;

  // X. Reprojection errors.
  std::unordered_map<Key, Eigen::Vector2d *> reproj_errs_;
  std::vector<Eigen::Vector2d> buf_reproj_errs_;
  std::vector<VectorReprojOnA *> reproj_errs_on_A_;
  std::vector<VectorReprojOnA> buf_reproj_errs_on_A_;
  std::vector<VectorReprojOnB *> reproj_errs_on_B_;
  std::vector<VectorReprojOnB> buf_reproj_errs_on_B_;

  // X. Linear system solver.
  std::vector<VectorForSchur *> vecs_for_schur_;
  std::vector<VectorForSchur> buf_vecs_schur_;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> schur_comp_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_, delta_A_, delta_B_;
};

} // namespace core

#endif