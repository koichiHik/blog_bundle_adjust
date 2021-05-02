
#ifndef __MATRIX_COMPUTATION_H__
#define __MATRIX_COMPUTATION_H__

// STL
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Eigen
#include <Eigen/Core>

// Original
#include <common_def.h>

namespace core {

// Camera Parameters, total 11
// tx, ty, tz, vx, vy, vz, fx, fy, cx, cy, s
//
// Matrix Aij : 2 x 11
void ComputeMatrixA(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    MatrixA &A);

// Point Parameters, total 3
// X, Y, Z
//
// Matrix Bij : 2 x 3
void ComputeMatrixB(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    MatrixB &B);

/**
 * Matrix V
 */
void ComputeMatrixW(const MatrixA &A, const MatrixB &B,
                    const Eigen::Matrix2d &Cov, MatrixW &W);

/**
 * Matrix U
 */
void ComputeMatrixU(const std::vector<MatrixA *> &As,
                    const std::vector<Eigen::Matrix2d *> &Covs, MatrixU &U);

/**
 * Matrix V
 */
void ComputeMatrixV(const std::vector<MatrixB *> &Bs,
                    const std::vector<Eigen::Matrix2d *> &Cov, MatrixV &V);

/**
 * Matrix Y
 */
void ComputeMatrixY(const std::unordered_map<Key, MatrixW *> &Ws,
                    const std::vector<MatrixV *> &Vs,
                    std::unordered_map<Key, MatrixY *> &Ys);

/**
 * Vector for Schur complement
 */
void ComputeVectorForSchur(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixY *> &Ys,
    const std::vector<VectorReprojOnA *> &reproj_error_on_As,
    const std::vector<VectorReprojOnB *> &reproj_error_on_Bs,
    std::vector<VectorForSchur *> &SchurVec);

/**
 * Schur complement.
 */
void ComputeSchurComplement(
    const std::vector<MatrixU *> &Us,
    const std::unordered_map<Key, MatrixY *> &Ys,
    const std::unordered_map<Key, MatrixW *> &Ws,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &SchurComp);

/**
 * Add dumping factor.
 */
void AddDumpingFactor(double myu, std::vector<MatrixU *> &Us,
                      std::vector<MatrixV *> &Vs);

/**
 * Delta computation.
 */
void ComputeDeltaA(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &SchurComp,
    const std::vector<VectorForSchur *> &SchurVec,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaA);

void ComputeDeltaB(const std::vector<MatrixV *> &Vs,
                   const std::vector<VectorReprojOnB *> &reproj_error_on_Bs,
                   const std::unordered_map<Key, MatrixW *> &Ws,
                   const Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaA,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> &deltaB);

/**
 * Reprojection error computation.
 */
void ComputeReprojectionError(
    const std::vector<Camera *> &cameras, const std::vector<Track *> &tracks,
    const std::vector<Eigen::Vector3d *> &points,
    std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error);

void ComputeReprojectionErrorOnCameraParameters(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixA *> &As,
    const std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error,
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixCov *>
        &conv_invs,
    std::vector<VectorReprojOnA *> &reproj_errs_on_params);

void ComputeReprojectionErrorOnPointParameters(
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixB *> &Bs,
    const std::unordered_map<std::pair<TrackId, CameraId>, Eigen::Vector2d *>
        &reproj_error,
    const std::unordered_map<std::pair<TrackId, CameraId>, MatrixCov *>
        &conv_invs,
    std::vector<VectorReprojOnB *> &reproj_errs_on_params);

} // namespace core

#endif