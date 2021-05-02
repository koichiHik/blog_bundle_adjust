

#ifndef __DERIVATIVE_V_SECOND_WITH_ANGLE_AXIS_H__
#define __DERIVATIVE_V_SECOND_WITH_ANGLE_AXIS_H__

// STL
#include <cmath>

// Eigen
#include <Eigen/Core>

// Original
#include "rotations.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

namespace optimization {

/***************************************************************/
/********************** Derivative V ***************************/
/***************************************************************/

inline double dv2_dXdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta,
                   2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           2 * K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) * pow(rot(1), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           2 * K(2, 2) * rot(1) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -2 * pow(K(2, 2), 2) * rot(0) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -2 * pow(K(2, 2), 2) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * rot(1) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dXdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
            rot(2) * sin(theta) / theta) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(2, 2) * rot(1) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dXdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dXdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
            rot(1) * sin(theta) / theta) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(2, 2) * rot(1) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dXds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dXdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dXdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
            rot(1) * sin(theta) / theta) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return K(1, 1) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dXdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 2) * K(2, 2) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-2 * K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -2 * pow(K(2, 2), 2) * p(1) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) -
                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 2) * K(2, 2) * p(0) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(1, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(0) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dXdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                rot(1) * rot(2) * sin(theta) / pow(theta, 3)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(1, 2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 1) * K(2, 2) * p(0) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(1, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dYdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -2 * pow(K(2, 2), 2) * rot(0) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2),
                   2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           2 * K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) * pow(rot(0), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * rot(0) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dYdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * rot(0) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dYdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dYdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * rot(0) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dYds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dYdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dYdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dYdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 2) * K(2, 2) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-2 * K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) +
                sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(1, 2) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) * K(2, 2) * p(0) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dYdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(1) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) * K(2, 2) * p(0) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dZdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -2 * pow(K(2, 2), 2) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * rot(1) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1,
                   2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           2 * K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dZdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dZdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dZdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dZds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dZdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dZdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dZdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       sin(theta) / theta) +
            K(1, 2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dZdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4)) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dfxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
            rot(2) * sin(theta) / theta) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(2, 2) * rot(1) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dfydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * rot(0) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dfydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dfydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dfydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 1.0 / (K(2, 2) * T(2) +
                  K(2, 2) * p(0) *
                      (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta) +
                  K(2, 2) * p(1) *
                      (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                  K(2, 2) * p(2) *
                      (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return 1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dfydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
           (T(1) +
            p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(2) * sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(2) * (-rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(2, 2) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dfydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                    rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(2, 2) * p(1) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           p(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dfydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4)) +
            p(2) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(2, 2) * p(0) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dfydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return p(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
            rot(1) * sin(theta) / theta) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(2, 2) * rot(1) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * rot(0) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dcydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           1.0 / (K(2, 2) * T(2) +
                  K(2, 2) * p(0) *
                      (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta) +
                  K(2, 2) * p(1) *
                      (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                  K(2, 2) * p(2) *
                      (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                    rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    sin(theta) / theta) +
            p(2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(2, 2) * p(1) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           p(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    sin(theta) / theta) +
            p(1) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(2, 2) * p(0) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           p(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dcydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(1) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 0;
  }
}

inline double dv2_dsdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dsdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
            rot(1) * sin(theta) / theta) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return K(1, 1) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (rot(0) * sin(theta) / theta +
            rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) *
           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 1.0 / (K(2, 2) * T(2) +
                  K(2, 2) * p(0) *
                      (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta) +
                  K(2, 2) * p(1) *
                      (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                  K(2, 2) * p(2) *
                      (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return 1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dtydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return 0;
  }
}

inline double dv2_dtzdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(1) * sin(theta) / theta)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 2) * K(2, 2) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dtzdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(1, 2) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 2) * K(2, 2) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dtzdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 2) * K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(1, 2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dtzdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtzdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
           (T(1) +
            p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(2) * sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(2) * (-rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(2, 2) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtzdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtzdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           1.0 / (K(2, 2) * T(2) +
                  K(2, 2) * p(0) *
                      (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta) +
                  K(2, 2) * p(1) *
                      (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                  K(2, 2) * p(2) *
                      (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(2, 2) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dtzds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtzdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dtzdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dtzdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -2 * K(1, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return -2 * K(1, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3);
  }
}

inline double dv2_dtzdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return -K(1, 2) * K(2, 2) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dtzdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return 2 * K(1, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3);
  }
}

inline double dv2_dtzdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return -K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -2 * pow(K(2, 2), 2) * p(1) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) +
                sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       sin(theta) / theta) +
            K(1, 2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                    rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(2, 2) * p(1) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           p(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dvxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                    rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    sin(theta) / theta) +
            p(2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(2, 2) * p(1) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           p(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dvxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return -K(1, 1) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return -K(1, 2) * K(2, 2) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(1, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(1, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(1, 2) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      3 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(2) * cos(theta) / pow(theta, 2) -
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * cos(theta) / pow(theta, 2) +
                      rot(1) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      5 * pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 3) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 3) * sin(theta) / pow(theta, 5) +
                      pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(0) * cos(theta) / pow(theta, 2) -
                      3 * rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) +
            p(2) *
                (K(1, 1) *
                     (pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 3) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 3) * sin(theta) / pow(theta, 5) +
                      pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      3 * rot(0) * cos(theta) / pow(theta, 2) +
                      3 * rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      5 * pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                 8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) -
                 3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                 3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                 6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(1) * cos(theta) / pow(theta, 2) +
                 rot(1) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) -
                 3 * pow(rot(0), 3) * cos(theta) / pow(theta, 4) +
                 3 * pow(rot(0), 3) * sin(theta) / pow(theta, 5) +
                 pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 3 * rot(0) * cos(theta) / pow(theta, 2) -
                 3 * rot(0) * sin(theta) / pow(theta, 3) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 4) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 4) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(0), 4) * (1 - cos(theta)) / pow(theta, 6) -
                 pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 5 * pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 10 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 2 * pow(K(2, 2), 2) * pow(p(1), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * p(1) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(1, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(1, 2) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3))) +
            p(2) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * cos(theta) / pow(theta, 2) +
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      8 * rot(0) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) -
                 3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) -
                 rot(0) * cos(theta) / pow(theta, 2) +
                 rot(0) * sin(theta) / pow(theta, 3) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) -
            K(2, 2) * p(1) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                 3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) +
                 3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                 rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(1) * cos(theta) / pow(theta, 2) -
                 rot(1) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                 5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                 8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 4 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 8 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 2) * K(2, 2) * p(0) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(1, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(1, 2) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(1, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4)))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                      8 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(2) * cos(theta) / pow(theta, 2) -
                      rot(2) * sin(theta) / pow(theta, 3))) +
            p(2) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(2) * cos(theta) / pow(theta, 2) +
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            2 * K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                 3 * rot(0) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) -
                 3 * rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 5) +
                 pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                 3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                 rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(2) * cos(theta) / pow(theta, 2) -
                 rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                     pow(theta, 5) -
                 8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 4 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 1) * K(2, 2) * p(0) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) -
                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                sin(theta) / theta) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       sin(theta) / theta)) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 2) * K(2, 2) * p(0) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(1, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(0) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(1, 2) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) * K(2, 2) * p(0) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(1, 1) + K(1, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(1, 1) * rot(0) + K(1, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4)) +
            p(2) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(2, 2) * p(0) * (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    sin(theta) / theta) +
            p(1) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(2, 2) * p(0) * (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           p(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dvyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return 2 * K(1, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3);
  }
}

inline double dv2_dvydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(1, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(1, 2) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3))) +
            p(2) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * cos(theta) / pow(theta, 2) +
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      8 * rot(0) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(1) *
                (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 sin(theta) / theta) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) -
                 3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) -
                 rot(0) * cos(theta) / pow(theta, 2) +
                 rot(0) * sin(theta) / pow(theta, 3) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) -
            K(2, 2) * p(1) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                 3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) +
                 3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                 rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(1) * cos(theta) / pow(theta, 2) -
                 rot(1) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                 5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                 8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 4 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 8 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 2) * K(2, 2) * p(0) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv2_dvydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(1, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(1, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(1, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      rot(2) * cos(theta) / pow(theta, 2) -
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(1), 3) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                      3 * rot(1) * cos(theta) / pow(theta, 2) +
                      3 * rot(1) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) +
            p(2) *
                (K(1, 1) *
                     (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(1), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      5 * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 sin(theta) / theta) -
            2 * K(2, 2) * p(1) *
                (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                 pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                 3 * pow(rot(1), 3) * cos(theta) / pow(theta, 4) -
                 3 * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                 3 * rot(1) * cos(theta) / pow(theta, 2) +
                 3 * rot(1) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                 3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) +
                 rot(0) * cos(theta) / pow(theta, 2) -
                 rot(0) * sin(theta) / pow(theta, 3) +
                 pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                 8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                 pow(rot(1), 4) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(1), 4) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(1), 4) * (1 - cos(theta)) / pow(theta, 6) -
                 5 * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 10 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                 2 * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -2 * K(1, 2) * K(2, 2) * pow(p(0), 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * pow(p(0), 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3);
  }
}

inline double dv2_dvydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(1, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                             pow(theta, 4) -
                         5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                             pow(theta, 5) +
                         8 * rot(0) * pow(rot(1), 2) * rot(2) *
                             (1 - cos(theta)) / pow(theta, 6) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         3 * rot(1) * pow(rot(2), 2) * cos(theta) /
                             pow(theta, 4) +
                         3 * rot(1) * pow(rot(2), 2) * sin(theta) /
                             pow(theta, 5) +
                         rot(1) * cos(theta) / pow(theta, 2) -
                         rot(1) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                             pow(theta, 4) -
                         5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                             pow(theta, 5) +
                         8 * rot(0) * rot(1) * pow(rot(2), 2) *
                             (1 - cos(theta)) / pow(theta, 6) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         3 * pow(rot(1), 2) * rot(2) * cos(theta) /
                             pow(theta, 4) -
                         3 * pow(rot(1), 2) * rot(2) * sin(theta) /
                             pow(theta, 5) -
                         rot(2) * cos(theta) / pow(theta, 2) +
                         rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2))) +
            p(2) *
                (K(1, 1) *
                     (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            2 * K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 3 * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                 3 * pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 5) -
                 rot(2) * cos(theta) / pow(theta, 2) +
                 rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (-rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 3 * rot(0) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) +
                 3 * rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 5) +
                 pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                     pow(theta, 4) +
                 5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                 4 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 1) * K(2, 2) * pow(p(0), 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                rot(1) * rot(2) * sin(theta) / pow(theta, 3)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(1, 2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return K(1, 1) * K(2, 2) * p(0) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(1, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dvzdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) *
               (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(1) * (1 - cos(theta)) / pow(theta, 2)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) * K(2, 2) * p(0) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) -
           K(2, 2) *
               (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4)) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (K(1, 1) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1));
  } else {
    return -K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvzdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (T(1) +
                p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                        rot(2) * sin(theta) / theta) +
                p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(2) * (-rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return p(0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv2_dvzdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvzdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(1) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (T(2) +
                p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                p(1) * (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 0;
  }
}

inline double dv2_dvzds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvzdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 0;
  } else {
    return 0;
  }
}

inline double dv2_dvzdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) *
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) /
           pow(K(2, 2) * T(2) +
                   K(2, 2) * p(0) *
                       (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(1) * sin(theta) / theta) +
                   K(2, 2) * p(1) *
                       (rot(0) * sin(theta) / theta +
                        rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                   K(2, 2) * p(2) *
                       (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1),
               2);
  } else {
    return 0;
  }
}

inline double dv2_dvzdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           K(2, 2) *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) -
           2 * K(2, 2) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3);
  } else {
    return -K(1, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(1, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(1, 2) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(1, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4)))) *
               (-K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                      8 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(2) * cos(theta) / pow(theta, 2) -
                      rot(2) * sin(theta) / pow(theta, 3))) +
            p(2) *
                (K(1, 1) *
                     (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(2) * cos(theta) / pow(theta, 2) +
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (-2 * K(2, 2) * p(0) *
                    (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(1) *
                    (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                     pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                     rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     sin(theta) / theta) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                     2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     2 * rot(0) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                 3 * rot(0) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) -
                 3 * rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 5) +
                 pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(1) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                 3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                 3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                 rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(2) * cos(theta) / pow(theta, 2) -
                 rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                     pow(theta, 5) -
                 8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 4 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return -K(1, 1) * K(2, 2) * p(0) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(1, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(1, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(1, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(1, 2) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) * (K(1, 1) *
                        (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                             pow(theta, 4) -
                         5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                             pow(theta, 5) +
                         8 * rot(0) * pow(rot(1), 2) * rot(2) *
                             (1 - cos(theta)) / pow(theta, 6) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         3 * rot(1) * pow(rot(2), 2) * cos(theta) /
                             pow(theta, 4) +
                         3 * rot(1) * pow(rot(2), 2) * sin(theta) /
                             pow(theta, 5) +
                         rot(1) * cos(theta) / pow(theta, 2) -
                         rot(1) * sin(theta) / pow(theta, 3)) +
                    K(1, 2) *
                        (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                             pow(theta, 4) -
                         5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                             pow(theta, 5) +
                         8 * rot(0) * rot(1) * pow(rot(2), 2) *
                             (1 - cos(theta)) / pow(theta, 6) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         3 * pow(rot(1), 2) * rot(2) * cos(theta) /
                             pow(theta, 4) -
                         3 * pow(rot(1), 2) * rot(2) * sin(theta) /
                             pow(theta, 5) -
                         rot(2) * cos(theta) / pow(theta, 2) +
                         rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2))) +
            p(2) *
                (K(1, 1) *
                     (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      4 * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (-2 * K(2, 2) * p(0) *
                    (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                     pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                     sin(theta) / theta) -
                2 * K(2, 2) * p(1) *
                    (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(2) * (1 - cos(theta)) / pow(theta, 2)) -
                2 * K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                     2 * rot(1) * (1 - cos(theta)) / pow(theta, 2))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 3 * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                 3 * pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 5) -
                 rot(2) * cos(theta) / pow(theta, 2) +
                 rot(2) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (-rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 3 * rot(0) * rot(1) * rot(2) * cos(theta) / pow(theta, 4) +
                 3 * rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 5) +
                 pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) +
                 8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                 (1 - cos(theta)) / pow(theta, 2)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                     pow(theta, 4) +
                 5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                 8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 2 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                 4 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return K(1, 1) * K(2, 2) * pow(p(0), 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double dv2_dvzdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(1, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(1, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(1, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * rot(2) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4)))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2) +
           (p(0) *
                (K(1, 1) *
                     (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 3) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                      3 * rot(2) * cos(theta) / pow(theta, 2) -
                      3 * rot(2) * sin(theta) / pow(theta, 3)) +
                 K(1, 2) *
                     (rot(0) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(0) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) -
                      rot(1) * cos(theta) / pow(theta, 2) +
                      rot(1) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(1, 1) *
                     (-pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(2), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(2), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      5 * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)) +
                 K(1, 2) *
                     (-rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4))) +
            p(2) *
                (K(1, 1) *
                     (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(1, 2) *
                     (-pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4)))) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) +
           (-2 * K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                 rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                 rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
            2 * K(2, 2) * p(1) *
                (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                 rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 4) +
                 rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
            2 * K(2, 2) * p(2) *
                (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4) -
                 pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                     pow(theta, 4))) *
               (-K(2, 2) * p(0) *
                    (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                     rot(1) * rot(2) * sin(theta) / pow(theta, 3)) -
                K(2, 2) * p(1) *
                    (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                     rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                     rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                     2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                         pow(theta, 4) +
                     rot(1) * (1 - cos(theta)) / pow(theta, 2)) -
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4) -
                     pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                     2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                         pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   3) +
           (-K(2, 2) * p(0) *
                (rot(0) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                 5 * rot(0) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                 8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                 6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                 rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                 3 * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                 3 * rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) -
                 rot(1) * cos(theta) / pow(theta, 2) +
                 rot(1) * sin(theta) / pow(theta, 3)) -
            K(2, 2) * p(1) *
                (-rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                 3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                 3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                 rot(0) * cos(theta) / pow(theta, 2) -
                 rot(0) * sin(theta) / pow(theta, 3) +
                 rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                 5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                 8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                     pow(theta, 6) +
                 3 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                 6 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) -
            K(2, 2) * p(2) *
                (-pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                 pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                 5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                     pow(theta, 5) -
                 8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                     pow(theta, 6) -
                 pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                 2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4))) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) /
                                       pow(theta, 2) +
                                   rot(2) * sin(theta) / theta) +
                        K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) /
                                       pow(theta, 2) -
                                   rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(1, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(1, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(1, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(1, 2) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1))) /
               pow(K(2, 2) * T(2) +
                       K(2, 2) * p(0) *
                           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                            rot(1) * sin(theta) / theta) +
                       K(2, 2) * p(1) *
                           (rot(0) * sin(theta) / theta + rot(1) * rot(2) *
                                                              (1 - cos(theta)) /
                                                              pow(theta, 2)) +
                       K(2, 2) * p(2) *
                           (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                            pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) +
                            1),
                   2);
  } else {
    return 0;
  }
}

} // namespace optimization

#pragma GCC diagnostic pop

#endif