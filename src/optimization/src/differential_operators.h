
#ifndef __DIFFERENTIAL_OPERATORS_H__
#define __DIFFERENTIAL_OPERATORS_H__

// STL
#include <map>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Original
#include "common_def.h"

namespace optimization {

void UpdateParameters(const Eigen::MatrixXd& delta_params,
                      const std::vector<Track>& tracks,
                      const std::map<size_t, size_t>& extrinsic_intrinsic_map,
                      std::vector<Eigen::Matrix3d>& K,
                      std::vector<Eigen::Vector3d>& T,
                      std::vector<Eigen::Vector3d>& Rot,
                      std::vector<Eigen::Vector3d>& points3d);

Eigen::MatrixXd ComputeResidualVector(
    const std::vector<Eigen::Matrix3d>& K,
    const std::vector<Eigen::Vector3d>& T,
    const std::vector<Eigen::Vector3d>& Rot,
    const std::vector<Eigen::Vector3d>& points3d,
    const std::vector<Track>& tracks,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map);

Eigen::MatrixXd CreateGaugeFixMatrix(size_t num_params, size_t num_constraints,
                                     const std::vector<size_t>& gauge_indices);

Eigen::MatrixXd ComputeGradient(
    const std::vector<Eigen::Matrix3d>& K,
    const std::vector<Eigen::Vector3d>& T,
    const std::vector<Eigen::Vector3d>& Rot,
    const std::vector<Eigen::Vector3d>& points3d,
    const std::vector<Track>& tracks,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map);

Eigen::MatrixXd ComputeJacobian(
    const std::vector<Eigen::Matrix3d>& K,
    const std::vector<Eigen::Vector3d>& T,
    const std::vector<Eigen::Vector3d>& Rot,
    const std::vector<Eigen::Vector3d>& points3d,
    const std::vector<Track>& tracks,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map);

Eigen::MatrixXd ComputeHessian(
    const std::vector<Eigen::Matrix3d>& K,
    const std::vector<Eigen::Vector3d>& T,
    const std::vector<Eigen::Vector3d>& Rot,
    const std::vector<Eigen::Vector3d>& points3d,
    const std::vector<Track>& tracks,
    const std::map<size_t, size_t>& extrinsic_intrinsic_map);

Eigen::MatrixXd MakeHessianPositiveDefiniteViaMultipleIdentity(
    const Eigen::MatrixXd& H, double multiple, double beta);

}  // namespace optimization

#endif