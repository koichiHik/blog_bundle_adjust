
#ifndef __GRADIENT_DESCENT_H__
#define __GRADIENT_DESCENT_H__

// X. STL
#include <vector>

// X. Eigen
#include <Eigen/Dense>

// X. Original
#include <common_def.h>

namespace core {

void LineSearchWithBackTracking(
    const std::vector<Track> &track_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map, size_t max_itr,
    double rho, double c, double &alpha_star);

void LineSearchWithStrongWolfe(
    const std::vector<Track> &tracks_src,
    const std::vector<Eigen::Matrix3d> &K_src,
    const std::vector<Eigen::Vector3d> &T_src,
    const std::vector<Eigen::Vector3d> &Rot_src,
    const std::vector<Eigen::Vector3d> &points3d_src,
    const std::map<size_t, size_t> &extrinsic_intrinsic_map, size_t max_itr,
    double alpha_max, double c1, double c2, double &alpha_star);

class GradientDescent {

public:
  void Optimize(const std::vector<Track> &tracks,
                const std::vector<Camera> &extrinsics,
                const std::vector<Eigen::Matrix3d> &intrinsics,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map,
                const std::vector<Eigen::Vector3d> &points3d,
                std::vector<Camera> &extrinsics_dst,
                std::vector<Eigen::Matrix3d> &intrinsics_dst,
                std::vector<Eigen::Vector3d> &points3d_dst);

private:
  const double INITIAL_STEP = 0.00000000000001;
  const double TERM_GRADIENT_NORM = 0.00000001;
  const double TERM_REPJ_DIFF = 1e-18;
  const double TERM_STEP_LENGTH = 1e-14;
};

class GradientDescentWithLineSearch {

public:
  void Optimize(const std::vector<Track> &tracks,
                const std::vector<Camera> &extrinsics,
                const std::vector<Eigen::Matrix3d> &intrinsics,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map,
                const std::vector<Eigen::Vector3d> &points3d,
                std::vector<Camera> &extrinsics_dst,
                std::vector<Eigen::Matrix3d> &intrinsics_dst,
                std::vector<Eigen::Vector3d> &points3d_dst);

private:
  const double INITIAL_STEP = 0.00000000000001;
  const double TERM_GRADIENT_NORM = 0.00000001;
  const double TERM_REPJ_DIFF = 1e-18;
  const double TERM_STEP_LENGTH = 1e-14;
};

class PureNewton {

public:
  void Optimize(const std::vector<Track> &tracks,
                const std::vector<Camera> &extrinsics,
                const std::vector<Eigen::Matrix3d> &intrinsics,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map,
                const std::vector<Eigen::Vector3d> &points3d,
                std::vector<Camera> &extrinsics_dst,
                std::vector<Eigen::Matrix3d> &intrinsics_dst,
                std::vector<Eigen::Vector3d> &points3d_dst);

private:
  const double TERM_DELTA_NORM = 0.00000001;
  const double TERM_REPJ_DIFF = 1e-10;
  const double BETA = 0.0001;
};

class CeresSolver {

public:
  void Optimize(const std::vector<Track> &tracks,
                const std::vector<Camera> &extrinsics,
                const std::vector<Eigen::Matrix3d> &intrinsics,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map,
                const std::vector<Eigen::Vector3d> &points3d,
                std::vector<Camera> &extrinsics_dst,
                std::vector<Eigen::Matrix3d> &intrinsics_dst,
                std::vector<Eigen::Vector3d> &points3d_dst);

private:
};

} // namespace core

#endif