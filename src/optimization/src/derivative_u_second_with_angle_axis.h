

#ifndef __DERIVATIVE_U_SECOND_WITH_ANGLE_AXIS_H__
#define __DERIVATIVE_U_SECOND_WITH_ANGLE_AXIS_H__

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
/********************** Derivative U ***************************/
/***************************************************************/

inline double du2_dXdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta,
                   2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
    return 2 * pow(K(2, 2), 2) * pow(rot(1), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           2 * K(2, 2) * rot(1) *
               (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
                   2) -
           K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
    return -2 * pow(K(2, 2), 2) * rot(0) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
                   2) -
           K(2, 2) *
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * rot(1) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
    return K(2, 2) * rot(1) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dXdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dXdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dXdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dXds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dXdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return K(0, 0) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dXdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return K(0, 1) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dXdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
    return K(0, 2) * K(2, 2) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(0, 1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) -
                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                sin(theta) / theta) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
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
    return -K(0, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(0) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dXdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                rot(1) * rot(2) * sin(theta) / pow(theta, 3)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(0, 2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
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
    return K(0, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           K(2, 2) * rot(1) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
                   2) -
           K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
    return -2 * pow(K(2, 2), 2) * rot(0) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2),
                   2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
    return 2 * pow(K(2, 2), 2) * pow(rot(0), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * rot(0) *
               (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
                   2) -
           K(2, 2) *
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
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
    return -K(2, 2) * rot(0) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dYdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dYdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dYdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dYds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dYdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return -K(0, 0) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dYdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return -K(0, 1) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dYdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
    return -K(0, 2) * K(2, 2) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) +
                sin(theta) / theta) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
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
    return K(0, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(0, 2) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
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
    return -2 * pow(K(2, 2), 2) * p(0) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dYdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(1) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
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
    return -K(0, 0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
           K(2, 2) * rot(0) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
                   2) -
           K(2, 2) *
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * rot(1) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
                   2) -
           K(2, 2) *
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * rot(0) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                       const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                       const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 * pow(K(2, 2), 2) *
               pow(-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1,
                   2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -K(2, 2) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dZdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dZdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dZdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dZds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dZdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return -K(0, 0) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dZdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return -K(0, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dZdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
    return -K(0, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       sin(theta) / theta) +
            K(0, 2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
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
    return -K(0, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(0, 1) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
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
    return K(0, 0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dZdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 1) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
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
    return -K(2, 2) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dfxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
    return K(2, 2) * rot(1) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           1.0 / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                  K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dfxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
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
    return -K(2, 2) * rot(0) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           rot(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dfxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) *
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
                   2) +
           (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -K(2, 2) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           rot(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                     K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dfxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
           (T(0) +
            p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(2) * sin(theta) / theta) +
            p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * sin(theta) / theta)) /
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
    return -K(2, 2) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dfxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4)) +
            p(1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -K(2, 2) * p(1) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dfxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return K(2, 2) * p(0) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           p(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dfxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    sin(theta) / theta) +
            p(2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -p(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                    K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dfydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dfydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dcydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dsdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return K(0, 0) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return -K(0, 0) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) *
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
    return -K(0, 0) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) /
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
    return -K(0, 0) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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
    return -K(0, 0) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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
    return K(0, 0) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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

inline double du2_dtydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return K(0, 1) * K(2, 2) * rot(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return -K(0, 1) * K(2, 2) * rot(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) *
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
    return -K(0, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) /
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
    return -K(0, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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
    return -K(0, 1) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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
    return K(0, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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

inline double du2_dtzdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) *
                    (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(2) * sin(theta) / theta) +
                K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
    return K(0, 2) * K(2, 2) * rot(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dtzdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                           rot(2) * sin(theta) / theta) +
                K(0, 1) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                K(0, 2) *
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
    return -K(0, 2) * K(2, 2) * rot(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dtzdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 2) * K(2, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                           rot(1) * sin(theta) / theta) +
                K(0, 1) * (-rot(0) * sin(theta) / theta +
                           rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(0, 2) *
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
    return -K(0, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dtzdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
           (T(0) +
            p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(2) * sin(theta) / theta) +
            p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * sin(theta) / theta)) /
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
    return -K(2, 2) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtzdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtzdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtzdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtzds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dtzdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 0) * K(2, 2) /
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
    return -K(0, 0) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtzdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(0, 1) * K(2, 2) /
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
    return -K(0, 1) * K(2, 2) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dtzdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -2 * K(0, 2) * K(2, 2) /
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -2 * K(0, 2) * K(2, 2) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3);
  }
}

inline double du2_dtzdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(0, 2) * K(2, 2) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dtzdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return K(0, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dtzdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(2, 2) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvxdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(0, 1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) +
                sin(theta) / theta) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
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
    return K(0, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                       rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       sin(theta) / theta) +
            K(0, 2) * (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
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
    return -K(0, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * p(1) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4)) +
            p(1) * (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                    rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -K(2, 2) * p(1) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvxdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvxdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvxdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvxds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvxdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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
    return -K(0, 0) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvxdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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
    return -K(0, 1) * K(2, 2) * p(1) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvxdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(0, 2) * K(2, 2) * p(1) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           2 * pow(K(2, 2), 2) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(0, 0) *
                         (-rot(0) * pow(rot(1), 2) * sin(theta) /
                              pow(theta, 3) +
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 1) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(1) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(0), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) * (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                sin(theta) / theta)) +
                p(2) *
                    (K(0, 0) *
                         (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                                rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                sin(theta) / theta) +
                     K(0, 2) *
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
                (K(0, 0) *
                     (-pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 3) * rot(1) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(0), 3) * rot(1) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(0), 3) * rot(1) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(0), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 5) +
                      3 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(2) * cos(theta) / pow(theta, 2) +
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 3) * rot(2) * cos(theta) / pow(theta, 4) -
                      5 * pow(rot(0), 3) * rot(2) * sin(theta) / pow(theta, 5) +
                      8 * pow(rot(0), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                      3 * pow(rot(0), 2) * rot(1) * cos(theta) / pow(theta, 4) +
                      3 * pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 5) +
                      3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           2 * K(2, 2) * p(1) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(0, 0) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) *
                (K(0, 0) *
                     (-rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -2 * pow(K(2, 2), 2) * p(0) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * p(1) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvxdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
                (K(0, 0) * (-rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) +
                            5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) -
                            8 * rot(0) * pow(rot(1), 2) * rot(2) *
                                (1 - cos(theta)) / pow(theta, 6) -
                            rot(0) * pow(rot(2), 3) * cos(theta) /
                                pow(theta, 4) +
                            5 * rot(0) * pow(rot(2), 3) * sin(theta) /
                                pow(theta, 5) -
                            8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                                pow(theta, 6) -
                            2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                            4 * rot(0) * rot(2) * (1 - cos(theta)) /
                                pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(2, 2) * p(1) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvydX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                    pow(theta, 4) -
                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                sin(theta) / theta) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
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
    return -K(0, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           2 * pow(K(2, 2), 2) * p(0) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) * rot(1) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           K(2, 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 1) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4)) +
            K(0, 2) * (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
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
    return -2 * pow(K(2, 2), 2) * p(0) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * rot(0) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) / pow(theta, 4) -
                pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(0, 1) * (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                       pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
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
    return K(0, 0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                    2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                    rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            p(2) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                    pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return K(2, 2) * p(0) * (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           p(2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dvydfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvydcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvydcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvyds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvydtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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
    return K(0, 0) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvydty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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
    return K(0, 1) * K(2, 2) * p(0) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvydtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return K(0, 2) * K(2, 2) * p(0) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           2 * pow(K(2, 2), 2) * p(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) -
           K(2, 2) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(0, 0) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) *
                (K(0, 0) *
                     (-rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      2 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                      4 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * pow(rot(1), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(1), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(1), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * pow(rot(1), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 5) +
                      rot(0) * cos(theta) / pow(theta, 2) -
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -2 * pow(K(2, 2), 2) * p(0) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           K(2, 2) * p(0) * (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) -
           K(2, 2) * p(1) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                          2 * pow(rot(1), 3) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 2) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * rot(1) * sin(theta) /
                              pow(theta, 3) +
                          2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                              pow(theta, 4) -
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4)) +
                     K(0, 2) *
                         (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                          pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                          2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
                (K(0, 0) *
                     (-pow(rot(1), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      5 * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (rot(0) * pow(rot(1), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(1), 2) * rot(2) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 5) -
                      rot(2) * cos(theta) / pow(theta, 2) +
                      rot(2) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) * (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) -
                            5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) +
                            8 * rot(0) * pow(rot(1), 2) * rot(2) *
                                (1 - cos(theta)) / pow(theta, 6) +
                            rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                            2 * rot(0) * rot(2) * (1 - cos(theta)) /
                                pow(theta, 4) -
                            pow(rot(1), 3) * sin(theta) / pow(theta, 3) -
                            3 * pow(rot(1), 3) * cos(theta) / pow(theta, 4) +
                            3 * pow(rot(1), 3) * sin(theta) / pow(theta, 5) +
                            3 * rot(1) * cos(theta) / pow(theta, 2) -
                            3 * rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return 2 * pow(K(2, 2), 2) * pow(p(0), 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   3) +
           2 * K(2, 2) * p(0) * (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvydvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(0, 0) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) *
                (K(0, 0) *
                     (-pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      8 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
                     (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 2) *
                     (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
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
                      rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(0, 0) * (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) -
                            5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) +
                            8 * rot(0) * pow(rot(1), 2) * rot(2) *
                                (1 - cos(theta)) / pow(theta, 6) +
                            rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                            2 * rot(0) * rot(2) * (1 - cos(theta)) /
                                pow(theta, 4) +
                            rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 3) +
                            3 * rot(1) * pow(rot(2), 2) * cos(theta) /
                                pow(theta, 4) -
                            3 * rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 5) -
                            rot(1) * cos(theta) / pow(theta, 2) +
                            rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) * (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                                pow(theta, 4) -
                            5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 5) +
                            8 * rot(0) * rot(1) * pow(rot(2), 2) *
                                (1 - cos(theta)) / pow(theta, 6) +
                            rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                            2 * rot(0) * rot(1) * (1 - cos(theta)) /
                                pow(theta, 4) -
                            pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 3) -
                            3 * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) +
                            3 * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) +
                            rot(2) * cos(theta) / pow(theta, 2) -
                            rot(2) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return K(2, 2) * p(0) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvzdX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * sin(theta) / theta) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                rot(1) * rot(2) * sin(theta) / pow(theta, 3)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                       sin(theta) / theta) +
            K(0, 2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
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
    return K(0, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           K(2, 2) * rot(1) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvzdY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (rot(0) * sin(theta) / theta +
                rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                rot(1) * (1 - cos(theta)) / pow(theta, 2)) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
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
           (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                       pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                           pow(theta, 4) -
                       pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                       2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                       2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
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
    return -K(0, 0) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
           K(2, 2) * rot(0) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double du2_dvzdZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                        const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                        const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
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
           (K(0, 0) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                       rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
            K(0, 1) * (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                       rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                       rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                       2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                           pow(theta, 4) +
                       rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
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
    return -K(2, 2) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvzdfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                    2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                    2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(1) * (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                        pow(theta, 4) -
                    pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                    pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    sin(theta) / theta) +
            p(2) * (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                    2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                        pow(theta, 4) +
                    rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
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
               (T(0) +
                p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                        pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
                p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                        rot(2) * sin(theta) / theta) +
                p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
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
    return -p(1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                    K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du2_dvzdfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvzdcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvzdcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvzds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du2_dvzdtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) *
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

inline double du2_dvzdty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) *
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

inline double du2_dvzdtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) *
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
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(2, 2) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvzdvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 1) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(1) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (pow(rot(0), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         sin(theta) / theta)) +
            p(2) * (K(0, 0) *
                        (pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(0), 2) * sin(theta) / pow(theta, 3) +
                         rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         sin(theta) / theta) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
                (K(0, 0) * (-rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) +
                            5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) -
                            8 * rot(0) * pow(rot(1), 2) * rot(2) *
                                (1 - cos(theta)) / pow(theta, 6) -
                            rot(0) * pow(rot(2), 3) * cos(theta) /
                                pow(theta, 4) +
                            5 * rot(0) * pow(rot(2), 3) * sin(theta) /
                                pow(theta, 5) -
                            8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                                pow(theta, 6) -
                            2 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                            4 * rot(0) * rot(2) * (1 - cos(theta)) /
                                pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * rot(1) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      3 * rot(0) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) -
                      3 * rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) -
                      rot(0) * cos(theta) / pow(theta, 2) +
                      rot(0) * sin(theta) / pow(theta, 3) +
                      rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (pow(rot(0), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * pow(rot(0), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * pow(rot(0), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      pow(rot(0), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                      3 * rot(0) * rot(1) * rot(2) * cos(theta) /
                          pow(theta, 4) +
                      3 * rot(0) * rot(1) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      2 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) +
                      (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return -K(2, 2) * p(1) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvzdvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
                        (-pow(rot(1), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta)) +
            p(1) * (K(0, 0) *
                        (rot(0) * pow(rot(1), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(1), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(1) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(1) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4)) +
                    K(0, 2) *
                        (rot(0) * rot(1) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(1), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-rot(0) * rot(1) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(1) * sin(theta) / pow(theta, 3) +
                         pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) * (K(0, 0) *
                        (-pow(rot(1), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 1) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                         sin(theta) / theta) +
                    K(0, 2) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) * (K(0, 0) *
                        (rot(0) * rot(1) * rot(2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * rot(1) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                         pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         sin(theta) / theta) +
                    K(0, 1) *
                        (-pow(rot(0), 2) * rot(2) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                             pow(theta, 4) -
                         pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                         2 * pow(rot(2), 3) * (1 - cos(theta)) / pow(theta, 4) -
                         2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
                        (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
            p(2) * (K(0, 0) *
                        (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                         rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                         rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                    K(0, 1) *
                        (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                         rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                         rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                         2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                             pow(theta, 4) +
                         rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                    K(0, 2) *
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
           (p(0) *
                (K(0, 0) *
                     (-pow(rot(1), 3) * rot(2) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(1), 3) * rot(2) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(1), 3) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      rot(1) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) +
                      5 * rot(1) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      8 * rot(1) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      4 * rot(1) * rot(2) * sin(theta) / pow(theta, 3) +
                      8 * rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 4)) +
                 K(0, 1) *
                     (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 2) *
                     (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
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
                      rot(2) * sin(theta) / pow(theta, 3))) +
            p(1) *
                (K(0, 0) * (rot(0) * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) -
                            5 * rot(0) * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) +
                            8 * rot(0) * pow(rot(1), 2) * rot(2) *
                                (1 - cos(theta)) / pow(theta, 6) +
                            rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                            2 * rot(0) * rot(2) * (1 - cos(theta)) /
                                pow(theta, 4) +
                            rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 3) +
                            3 * rot(1) * pow(rot(2), 2) * cos(theta) /
                                pow(theta, 4) -
                            3 * rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 5) -
                            rot(1) * cos(theta) / pow(theta, 2) +
                            rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) * (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                                pow(theta, 4) -
                            5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                                pow(theta, 5) +
                            8 * rot(0) * rot(1) * pow(rot(2), 2) *
                                (1 - cos(theta)) / pow(theta, 6) +
                            rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                            2 * rot(0) * rot(1) * (1 - cos(theta)) /
                                pow(theta, 4) -
                            pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 3) -
                            3 * pow(rot(1), 2) * rot(2) * cos(theta) /
                                pow(theta, 4) +
                            3 * pow(rot(1), 2) * rot(2) * sin(theta) /
                                pow(theta, 5) +
                            rot(2) * cos(theta) / pow(theta, 2) -
                            rot(2) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
    return K(2, 2) * p(0) * (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                   K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
               2);
  }
}

inline double du2_dvzdvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                         const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                         const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return 2 *
               (p(0) *
                    (K(0, 0) * (-pow(rot(1), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(1), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 1) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) +
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) -
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                                sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3))) +
                p(1) *
                    (K(0, 0) * (rot(0) * rot(1) * rot(2) * sin(theta) /
                                    pow(theta, 3) -
                                2 * rot(0) * rot(1) * rot(2) *
                                    (1 - cos(theta)) / pow(theta, 4) -
                                pow(rot(2), 2) * cos(theta) / pow(theta, 2) +
                                pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                                sin(theta) / theta) +
                     K(0, 1) * (-pow(rot(0), 2) * rot(2) * sin(theta) /
                                    pow(theta, 3) +
                                2 * pow(rot(0), 2) * rot(2) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                                2 * pow(rot(2), 3) * (1 - cos(theta)) /
                                    pow(theta, 4) -
                                2 * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
                         (rot(0) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(0) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(0) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * rot(2) * cos(theta) / pow(theta, 2) -
                          rot(1) * rot(2) * sin(theta) / pow(theta, 3)) +
                     K(0, 1) *
                         (-rot(0) * rot(2) * cos(theta) / pow(theta, 2) +
                          rot(0) * rot(2) * sin(theta) / pow(theta, 3) +
                          rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                          2 * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                              pow(theta, 4) +
                          rot(1) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
                (K(0, 0) *
                     (-pow(rot(1), 2) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) +
                      5 * pow(rot(1), 2) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) -
                      8 * pow(rot(1), 2) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) -
                      pow(rot(1), 2) * sin(theta) / pow(theta, 3) +
                      2 * pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      pow(rot(2), 4) * cos(theta) / pow(theta, 4) +
                      5 * pow(rot(2), 4) * sin(theta) / pow(theta, 5) -
                      8 * pow(rot(2), 4) * (1 - cos(theta)) / pow(theta, 6) -
                      5 * pow(rot(2), 2) * sin(theta) / pow(theta, 3) +
                      10 * pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 4) -
                      2 * (1 - cos(theta)) / pow(theta, 2)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (rot(0) * rot(1) * pow(rot(2), 2) * cos(theta) /
                          pow(theta, 4) -
                      5 * rot(0) * rot(1) * pow(rot(2), 2) * sin(theta) /
                          pow(theta, 5) +
                      8 * rot(0) * rot(1) * pow(rot(2), 2) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      rot(0) * rot(1) * sin(theta) / pow(theta, 3) -
                      2 * rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 4) +
                      pow(rot(2), 3) * sin(theta) / pow(theta, 3) +
                      3 * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                      3 * pow(rot(2), 3) * sin(theta) / pow(theta, 5) -
                      3 * rot(2) * cos(theta) / pow(theta, 2) +
                      3 * rot(2) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
                (K(0, 0) *
                     (rot(0) * pow(rot(2), 3) * cos(theta) / pow(theta, 4) -
                      5 * rot(0) * pow(rot(2), 3) * sin(theta) / pow(theta, 5) +
                      8 * rot(0) * pow(rot(2), 3) * (1 - cos(theta)) /
                          pow(theta, 6) +
                      3 * rot(0) * rot(2) * sin(theta) / pow(theta, 3) -
                      6 * rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 4) -
                      rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 3) -
                      3 * rot(1) * pow(rot(2), 2) * cos(theta) / pow(theta, 4) +
                      3 * rot(1) * pow(rot(2), 2) * sin(theta) / pow(theta, 5) +
                      rot(1) * cos(theta) / pow(theta, 2) -
                      rot(1) * sin(theta) / pow(theta, 3)) +
                 K(0, 1) *
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
                 K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) *
                    (K(0, 0) *
                         (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 1) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(2) * sin(theta) / theta) +
                     K(0, 2) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(1) * sin(theta) / theta)) +
                p(1) *
                    (K(0, 0) *
                         (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                          rot(2) * sin(theta) / theta) +
                     K(0, 1) *
                         (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                          pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) +
                          1) +
                     K(0, 2) *
                         (rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) +
                p(2) *
                    (K(0, 0) *
                         (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                          rot(1) * sin(theta) / theta) +
                     K(0, 1) *
                         (-rot(0) * sin(theta) / theta +
                          rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                     K(0, 2) *
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