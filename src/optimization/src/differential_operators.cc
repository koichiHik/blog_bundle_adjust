
// STL
#include <map>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Google
#include <glog/logging.h>

// Original
#include "common_def.h"
#include "derivative_first_with_angle_axis.h"
#include "derivative_u_second_with_angle_axis.h"
#include "derivative_v_second_with_angle_axis.h"
#include "differential_operators.h"

namespace optimization {

void UpdateParameters(const Eigen::MatrixXd &delta_params,
                      const std::vector<Track> &tracks,
                      const std::map<size_t, size_t> &extrinsic_intrinsic_map,
                      std::vector<Eigen::Matrix3d> &K,
                      std::vector<Eigen::Vector3d> &T,
                      std::vector<Eigen::Vector3d> &Rot,
                      std::vector<Eigen::Vector3d> &points3d) {

  // X. Update camera extrinsics.
  size_t cam_ext_num = T.size();
  size_t cam_ext_var_num = 6;
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {
    size_t idx = cam_ext_idx * cam_ext_var_num;

    // X. Camera translation.
    T[cam_ext_idx](0) += delta_params(idx, 0);
    T[cam_ext_idx](1) += delta_params(idx + 1, 0);
    T[cam_ext_idx](2) += delta_params(idx + 2, 0);

    // X. Camera rotation.
    Rot[cam_ext_idx](0) += delta_params(idx + 3, 0);
    Rot[cam_ext_idx](1) += delta_params(idx + 4, 0);
    Rot[cam_ext_idx](2) += delta_params(idx + 5, 0);
  }

  // X. Update camera intrinsics.
  size_t cam_int_num = K.size();
  size_t cam_int_var_num = 5;
  for (size_t cam_int_idx = 0; cam_int_idx < cam_int_num; cam_int_idx++) {
    size_t idx = cam_ext_num * cam_ext_var_num + cam_int_idx * cam_int_var_num;

    // X. Camera intrinsic.
    K[cam_int_idx](0, 0) += delta_params(idx, 0);
    K[cam_int_idx](1, 1) += delta_params(idx + 1, 0);
    K[cam_int_idx](0, 2) += delta_params(idx + 2, 0);
    K[cam_int_idx](1, 2) += delta_params(idx + 3, 0);
    K[cam_int_idx](0, 1) += delta_params(idx + 4, 0);
  }

  // X. Point locations.
  size_t trk_num = tracks.size();
  size_t trk_var_num = 3;
  for (size_t trk_idx = 0; trk_idx < trk_num; trk_idx++) {
    size_t idx = trk_var_num * trk_idx + cam_ext_num * cam_ext_var_num +
                 cam_int_num * cam_int_var_num;
    points3d[trk_idx](0) += delta_params(idx, 0);
    points3d[trk_idx](1) += delta_params(idx + 1, 0);
    points3d[trk_idx](2) += delta_params(idx + 2, 0);
  }
}

Eigen::MatrixXd
ComputeResidualVector(const std::vector<Eigen::Matrix3d> &K,
                      const std::vector<Eigen::Vector3d> &T,
                      const std::vector<Eigen::Vector3d> &Rot,
                      const std::vector<Eigen::Vector3d> &points3d,
                      const std::vector<Track> &tracks,
                      const std::map<size_t, size_t> &extrinsic_intrinsic_map) {

  // X. Prepare buffer.
  size_t cam_ext_num = T.size();
  size_t track_num = tracks.size();
  size_t projected_pnt_num = 0;
  for (const Track &track : tracks) {
    projected_pnt_num += track.size();
  }
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(projected_pnt_num * 2, 1);

  // X. Compose camera matrix.
  std::vector<Camera> cameras;
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {

    // X. Create rotation matrix.
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);

    // X. Create camera matrix.
    Camera cam;
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T[cam_ext_idx];
    cam = K[extrinsic_intrinsic_map.at(cam_ext_idx)] * cam;
    cameras.push_back(cam);
  }

  // X.
  size_t row_idx = 0;
  for (size_t trk_idx = 0; trk_idx < track_num; trk_idx++) {
    for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {
      // X. Only visible track.
      if (tracks[trk_idx].count(cam_ext_idx)) {
        Eigen::Vector2d meas = tracks[trk_idx].at(cam_ext_idx);
        Eigen::Vector2d proj =
            (cameras[cam_ext_idx] * points3d[trk_idx].homogeneous())
                .hnormalized();
        res.block<2, 1>(row_idx, 0) = proj - meas;
        row_idx += 2;
      }
    }
  }
  CHECK(row_idx == 2 * projected_pnt_num)
      << "Final row_idx is not valid. row_idx : " << row_idx
      << ", projected_pnt_num : " << projected_pnt_num;
  return res;
}

Eigen::MatrixXd
ComputeGradient(const std::vector<Eigen::Matrix3d> &K,
                const std::vector<Eigen::Vector3d> &T,
                const std::vector<Eigen::Vector3d> &Rot,
                const std::vector<Eigen::Vector3d> &points3d,
                const std::vector<Track> &tracks,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map) {

  CHECK(points3d.size() == tracks.size())
      << "Number of points and tracks must be the same.";

  // X. Create Gradient.
  size_t cam_ext_num = T.size();
  size_t cam_int_num = K.size();
  size_t track_num = tracks.size();

  size_t cam_ext_var_num = 6;
  size_t cam_int_var_num = 5;
  size_t points_var_num = 3;
  size_t len = cam_ext_var_num * cam_ext_num + cam_int_var_num * cam_int_num +
               points_var_num * track_num;
  Eigen::MatrixXd grad(len, 1);
  grad.setZero();

  std::vector<Camera> cameras;
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {

    // X. Create rotation matrix.
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);

    // X. Create camera matrix.
    Camera cam;
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T[cam_ext_idx];
    cam = K[extrinsic_intrinsic_map.at(cam_ext_idx)] * cam;
    cameras.push_back(cam);
  }

  // X. Gradient for camera extrinsic part.
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {

    for (size_t trk_idx = 0; trk_idx < track_num; trk_idx++) {

      // X. Only visible track.
      if (tracks[trk_idx].count(cam_ext_idx)) {

        Eigen::Vector2d meas = tracks[trk_idx].at(cam_ext_idx);
        Eigen::Vector3d x =
            cameras[cam_ext_idx] * points3d[trk_idx].homogeneous();
        Eigen::Vector2d proj = x.hnormalized();

        double u_res = proj(0) - meas(0);
        double v_res = proj(1) - meas(1);

        size_t cam_int_idx = extrinsic_intrinsic_map.at(cam_ext_idx);
        {
#if 1
          size_t col_idx = cam_ext_idx * cam_ext_var_num;
          grad(col_idx, 0) +=
              u_res * du_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

          grad(col_idx + 1, 0) +=
              u_res * du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

          grad(col_idx + 2, 0) +=
              u_res * du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

          grad(col_idx + 3, 0) +=
              u_res * du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

          grad(col_idx + 4, 0) +=
              u_res * du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

          grad(col_idx + 5, 0) +=
              u_res * du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x) +
              v_res * dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                             points3d[trk_idx], x);

#endif
        }

        // X. Gradient for camera intrinsic part.
        {
          double scale = 1.0;
          size_t col_idx =
              cam_ext_num * cam_ext_var_num + cam_int_idx * cam_int_var_num;
          grad(col_idx, 0) +=
              scale * (u_res * du_dfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x) +
                       v_res * dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x));

          grad(col_idx + 1, 0) +=
              scale * (u_res * du_dfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x) +
                       v_res * dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x));

          grad(col_idx + 2, 0) +=
              scale * (u_res * du_dcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x) +
                       v_res * dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x));

          grad(col_idx + 3, 0) +=
              scale * (u_res * du_dcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x) +
                       v_res * dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x));

          grad(col_idx + 4, 0) +=
              scale * (u_res * du_ds(T[cam_ext_idx], Rot[cam_ext_idx],
                                     K[cam_int_idx], points3d[trk_idx], x) +
                       v_res * dv_ds(T[cam_ext_idx], Rot[cam_ext_idx],
                                     K[cam_int_idx], points3d[trk_idx], x));
        }
        // X. Gradient for point part.
        {
#if 1
          size_t col_idx = cam_ext_num * cam_ext_var_num +
                           cam_int_num * cam_int_var_num +
                           trk_idx * points_var_num;

          grad(col_idx, 0) +=
              u_res * du_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x) +
              v_res * dv_dX(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x);

          grad(col_idx + 1, 0) +=
              u_res * du_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x) +
              v_res * dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x);

          grad(col_idx + 2, 0) +=
              u_res * du_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x) +
              v_res * dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                            points3d[trk_idx], x);
#endif
        }
      }
    }
  }

  {
    size_t cam_ext_idx = 0;
    grad(cam_ext_idx, 0) = 0.0;
    grad(cam_ext_idx + 1, 0) = 0.0;
    grad(cam_ext_idx + 2, 0) = 0.0;
    grad(cam_ext_idx + 3, 0) = 0.0;
    grad(cam_ext_idx + 4, 0) = 0.0;
    grad(cam_ext_idx + 5, 0) = 0.0;

    size_t cam_2nd_idx = cam_ext_var_num;
    grad(cam_2nd_idx, 0) = 0.0;
  }

  return grad;
}

Eigen::MatrixXd
ComputeJacobian(const std::vector<Eigen::Matrix3d> &K,
                const std::vector<Eigen::Vector3d> &T,
                const std::vector<Eigen::Vector3d> &Rot,
                const std::vector<Eigen::Vector3d> &points3d,
                const std::vector<Track> &tracks,
                const std::map<size_t, size_t> &extrinsic_intrinsic_map) {

  CHECK(points3d.size() == tracks.size())
      << "Number of points and tracks must be the same.";

  // X. Element count.
  size_t cam_ext_num = T.size();
  size_t cam_int_num = K.size();
  size_t track_num = tracks.size();
  size_t projected_pnt_num = 0;
  for (const Track &track : tracks) {
    projected_pnt_num += track.size();
  }

  // X. Variable count.
  size_t cam_ext_var_num = 6;
  size_t cam_int_var_num = 5;
  size_t points_var_num = 3;

  // X. Prepare buffer.
  size_t j_row = projected_pnt_num * 2;
  size_t j_col = cam_ext_num * cam_ext_var_num + cam_int_num * cam_int_var_num +
                 track_num * points_var_num;
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(j_row, j_col);

  std::vector<Camera> cameras;
  for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {

    // X. Create rotation matrix.
    Eigen::Matrix3d R;
    ConvertAngleAxisToRotationMatrix(Rot[cam_ext_idx], R);

    // X. Create camera matrix.
    Camera cam;
    cam.block<3, 3>(0, 0) = R;
    cam.block<3, 1>(0, 3) = T[cam_ext_idx];
    cam = K[extrinsic_intrinsic_map.at(cam_ext_idx)] * cam;
    cameras.push_back(cam);
  }

  // X. Compute jacobian.
  size_t row_idx = 0;
  for (size_t trk_idx = 0; trk_idx < track_num; trk_idx++) {

    for (size_t cam_ext_idx = 0; cam_ext_idx < cam_ext_num; cam_ext_idx++) {

      // X. Only visible track.
      if (tracks[trk_idx].count(cam_ext_idx)) {

        Eigen::Vector3d x =
            cameras[cam_ext_idx] * points3d[trk_idx].homogeneous();

        // X. Compute extrinsic dependency.
        {
          size_t col_idx = cam_ext_idx * cam_ext_var_num;
          size_t cam_int_idx = extrinsic_intrinsic_map.at(cam_ext_idx);

          // X. tx
          J(row_idx, col_idx) = du_dtx(T[cam_ext_idx], Rot[cam_ext_idx],
                                       K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx) =
              dv_dtx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. ty
          J(row_idx, col_idx + 1) =
              du_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 1) =
              dv_dty(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. tz
          J(row_idx, col_idx + 2) =
              du_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 2) =
              dv_dtz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          // X. vx
          J(row_idx, col_idx + 3) =
              du_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 3) =
              dv_dvx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. vy
          J(row_idx, col_idx + 4) =
              du_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 4) =
              dv_dvy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. vz
          J(row_idx, col_idx + 5) =
              du_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 5) =
              dv_dvz(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
        }

        // X. Compute intrinsic dependency
        {
          size_t cam_int_idx = extrinsic_intrinsic_map.at(cam_ext_idx);
          size_t col_idx =
              cam_ext_num * cam_ext_var_num + cam_int_idx * cam_int_var_num;

          // X. fx
          J(row_idx, col_idx) = du_dfx(T[cam_ext_idx], Rot[cam_ext_idx],
                                       K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx) =
              dv_dfx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. fy
          J(row_idx, col_idx + 1) =
              du_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 1) =
              dv_dfy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. cx
          J(row_idx, col_idx + 2) =
              du_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 2) =
              dv_dcx(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. cy
          J(row_idx, col_idx + 3) =
              du_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 3) =
              dv_dcy(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                     points3d[trk_idx], x);

          // X. s
          J(row_idx, col_idx + 4) = du_ds(T[cam_ext_idx], Rot[cam_ext_idx],
                                          K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 4) =
              dv_ds(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                    points3d[trk_idx], x);
        }

        // X. Compute point coords dependency.
        {
          size_t col_idx = cam_ext_num * cam_ext_var_num +
                           cam_int_num * cam_int_var_num +
                           trk_idx * points_var_num;

          size_t cam_int_idx = extrinsic_intrinsic_map.at(cam_ext_idx);

          // X.
          J(row_idx, col_idx) = du_dX(T[cam_ext_idx], Rot[cam_ext_idx],
                                      K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx) = dv_dX(T[cam_ext_idx], Rot[cam_ext_idx],
                                          K[cam_int_idx], points3d[trk_idx], x);

          // Y.
          J(row_idx, col_idx + 1) = du_dY(T[cam_ext_idx], Rot[cam_ext_idx],
                                          K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 1) =
              dv_dY(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                    points3d[trk_idx], x);

          // Z.
          J(row_idx, col_idx + 2) = du_dZ(T[cam_ext_idx], Rot[cam_ext_idx],
                                          K[cam_int_idx], points3d[trk_idx], x);
          J(row_idx + 1, col_idx + 2) =
              dv_dZ(T[cam_ext_idx], Rot[cam_ext_idx], K[cam_int_idx],
                    points3d[trk_idx], x);
        }

        row_idx += 2;
      }
    }
  }

  CHECK(row_idx == 2 * projected_pnt_num)
      << "Final row_idx is not valid. row_idx : " << row_idx
      << ", projected_pnt_num : " << projected_pnt_num;

  return J;
}

} // namespace optimization