

// STL

// GLOG
#include <glog/logging.h>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Google Test
#include <gtest/gtest.h>

// Original
#include "rotations.h"
#include <geometry.h>

TEST(RQ_DECOMP, NORMAL) {
  LOG(INFO) << "RQ Decomposition Test.";

  Eigen::Matrix<double, 3, 4> M;
#if 0
  M << 3.53553e2, 3.39645e2, 2.77744e2, -1.44946e6, -1.03528e2, 2.33212e1,
      4.59607e2, -6.32525e5, 7.07107e-1, -3.53553e-1, 6.12372e-1, -9.18559e2;
#endif
  M << 3.9923568756416135, 39.41768098301378, -0.7632898797149192,
      3.9591755089132286, -14.430231011327074, -0.9414415802377172,
      -27.450970108566686, -14.429433437768129, 0.012249240354938502,
      -0.00014574603756147602, -0.0005693070873097415, 0.012249358697517865;

  Eigen::Matrix3d K, R;
  core::ComputeInternalCalibration(M, K, R);

  double residue = (M.block<3, 3>(0, 0) - K * R).norm();
  CHECK(residue < 1e-10);

  LOG(INFO) << "M";
  LOG(INFO) << M(0, 0) << ", " << M(0, 1) << ", " << M(0, 2) << ", " << M(0, 3);
  LOG(INFO) << M(1, 0) << ", " << M(1, 1) << ", " << M(1, 2) << ", " << M(1, 3);
  LOG(INFO) << M(2, 0) << ", " << M(2, 1) << ", " << M(2, 2) << ", " << M(2, 3);

  LOG(INFO) << "K";
  LOG(INFO) << K(0, 0) << ", " << K(0, 1) << ", " << K(0, 2);
  LOG(INFO) << K(1, 0) << ", " << K(1, 1) << ", " << K(1, 2);
  LOG(INFO) << K(2, 0) << ", " << K(2, 1) << ", " << K(2, 2);

  LOG(INFO) << "Rot";
  LOG(INFO) << R(0, 0) << ", " << R(0, 1) << ", " << R(0, 2);
  LOG(INFO) << R(1, 0) << ", " << R(1, 1) << ", " << R(1, 2);
  LOG(INFO) << R(2, 0) << ", " << R(2, 1) << ", " << R(2, 2);

  LOG(INFO) << "R.determinant() : " << R.determinant();
  LOG(INFO) << "Residue : " << residue;
}

TEST(ANGLE_AXIS, NORMAL) {

  double rx = 30.0 / 180.0 * M_PI;
  double ry = 45.0 / 180.0 * M_PI;
  double rz = 60.0 / 180.0 * M_PI;

  Eigen::Matrix3d R_old, R_new, R_calc, R_via_angle_axis;
  {
    R_old = Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
            Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ());

    Eigen::AngleAxisd rot_axis(R_old);
    double theta = rot_axis.angle();
    double vx = rot_axis.axis()(0), vy = rot_axis.axis()(1),
           vz = rot_axis.axis()(2);
    R_new = rot_axis.toRotationMatrix();

    Eigen::Matrix3d W;
    W << 0, -vz, vy, vz, 0, -vx, -vy, vx, 0;
    R_calc =
        Eigen::Matrix3d::Identity() + sin(theta) * W + (1 - cos(theta)) * W * W;

    Eigen::AngleAxisd rot_axis_via_vec(theta, Eigen::Vector3d(vx, vy, vz));
    R_via_angle_axis = rot_axis_via_vec.toRotationMatrix();

    rot_axis_via_vec.axis().norm();
  }

  CHECK((R_new - R_old).norm() < DBL_EPSILON * 10.0);
  CHECK((R_calc - R_old).norm() < DBL_EPSILON * 10.0);
  CHECK((R_new - R_via_angle_axis).norm() < DBL_EPSILON * 10.0);
}

TEST(ANGLE_AXIS, CONVERT) {

  double rx = 30.0 / 180.0 * M_PI;
  double ry = 45.0 / 180.0 * M_PI;
  double rz = 60.0 / 180.0 * M_PI;

  Eigen::Matrix3d R_in, R_out;
  R_in = Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()) *
         Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()) *
         Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ());

  Eigen::Vector3d angle_axis;
  core::ConvertRotationMatrixToAngleAxis(R_in, angle_axis);
  core::ConvertAngleAxisToRotationMatrix(angle_axis, R_out);

  CHECK((R_in - R_out).norm() < DBL_EPSILON * 100);
}

TEST(CALIB_DECOMP, TEST) {

  core::Camera cam;
  cam << 3.9923568756416135, 39.41768098301378, -0.7632898797149192,
      3.9591755089132286, -14.430231011327074, -0.9414415802377172,
      -27.450970108566686, -14.429433437768129, 0.012249240354938502,
      -0.00014574603756147602, -0.0005693070873097415, 0.012249358697517865;

  Eigen::Matrix3d K, R_in, R_out;
  core::ComputeInternalCalibration(cam, K, R_in);

  // R_in.col(2) = -1 * R_in.col(2);
  LOG(INFO) << R_in.determinant();

  Eigen::Vector3d angle_axis;
  core::ConvertRotationMatrixToAngleAxis(R_in, angle_axis);
  core::ConvertAngleAxisToRotationMatrix(angle_axis, R_out);

  /*
  std::cout << R_in << std::endl;
  std::cout << R_out << std::endl;
  std::cout << (R_in - R_out).norm() << std::endl;
  */
}