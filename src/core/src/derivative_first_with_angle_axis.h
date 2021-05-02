
// STL
#include <cmath>

// Eigen
#include <Eigen/Core>

// Original
#include "derivative_common.h"
#include "rotations.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

namespace optimization {

/***************************************************************/
/********************** Derivative U ***************************/
/***************************************************************/

inline double du_dX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(0, 0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(0, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) /
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
    return K(2, 2) * rot(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(0, 0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(2) * sin(theta) / theta) +
            K(0, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(0, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
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
    return -K(2, 2) * rot(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(0, 0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(1) * sin(theta) / theta) +
            K(0, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(0, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
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
    return -K(2, 2) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (T(0) +
            p(0) * (-pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(2) * sin(theta) / theta) +
            p(2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(1) * sin(theta) / theta)) /
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
    return (T(0) + p(0) - p(1) * rot(2) + p(2) * rot(1)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

inline double du_dfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du_dcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (T(2) +
            p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * sin(theta) / theta) +
            p(1) * (rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
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
    return (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

inline double du_dcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double du_ds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (T(1) +
            p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(2) * sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(2) * (-rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
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
    return (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

inline double du_dtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 0) /
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
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 1) /
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
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(0, 2) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) -
           K(2, 2) *
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
    return K(0, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
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

inline double du_dvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
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
    return -K(2, 2) * p(1) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (-K(0, 1) * p(2) + K(0, 2) * p(1)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
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
    return K(2, 2) * p(0) *
               (K(0, 0) * T(0) + K(0, 1) * T(1) + K(0, 2) * T(2) +
                p(0) * (K(0, 0) + K(0, 1) * rot(2) - K(0, 2) * rot(1)) +
                p(1) * (-K(0, 0) * rot(2) + K(0, 1) + K(0, 2) * rot(0)) +
                p(2) * (K(0, 0) * rot(1) - K(0, 1) * rot(0) + K(0, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (K(0, 0) * p(2) - K(0, 2) * p(0)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double du_dvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(0, 0) *
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
    return (-K(0, 0) * p(1) + K(0, 1) * p(0)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

/***************************************************************/
/********************** Derivative V ***************************/
/***************************************************************/

inline double dv_dX(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(1, 1) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                       rot(2) * sin(theta) / theta) +
            K(1, 2) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                       rot(1) * sin(theta) / theta)) /
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
    return K(2, 2) * rot(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (K(1, 1) * rot(2) - K(1, 2) * rot(1)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv_dY(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(1, 1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            K(1, 2) * (rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
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
    return -K(2, 2) * rot(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (K(1, 1) + K(1, 2) * rot(0)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv_dZ(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                    const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                    const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return -K(2, 2) *
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
                   2) +
           (K(1, 1) * (-rot(0) * sin(theta) / theta +
                       rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            K(1, 2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                       pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
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
    return -K(2, 2) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (-K(1, 1) * rot(0) + K(1, 2)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv_dfx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double dv_dfy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (T(1) +
            p(0) * (rot(0) * rot(1) * (1 - cos(theta)) / pow(theta, 2) +
                    rot(2) * sin(theta) / theta) +
            p(1) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(2), 2) * (1 - cos(theta)) / pow(theta, 2) + 1) +
            p(2) * (-rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2))) /
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
    return (T(1) + p(0) * rot(2) + p(1) - p(2) * rot(0)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

inline double dv_dcx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double dv_dcy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (T(2) +
            p(0) * (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                    rot(1) * sin(theta) / theta) +
            p(1) * (rot(0) * sin(theta) / theta +
                    rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
            p(2) * (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                    pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) /
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
    return (T(2) - p(0) * rot(1) + p(1) * rot(0) + p(2)) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

inline double dv_ds(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double dv_dtx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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

inline double dv_dty(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 1) /
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
    return K(1, 1) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv_dtz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return K(1, 2) /
               (K(2, 2) * T(2) +
                K(2, 2) * p(0) *
                    (rot(0) * rot(2) * (1 - cos(theta)) / pow(theta, 2) -
                     rot(1) * sin(theta) / theta) +
                K(2, 2) * p(1) *
                    (rot(0) * sin(theta) / theta +
                     rot(1) * rot(2) * (1 - cos(theta)) / pow(theta, 2)) +
                K(2, 2) * p(2) *
                    (-pow(rot(0), 2) * (1 - cos(theta)) / pow(theta, 2) -
                     pow(rot(1), 2) * (1 - cos(theta)) / pow(theta, 2) + 1)) -
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
    return K(1, 2) / (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                      K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) -
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

inline double dv_dvx(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
                         2 * rot(0) * (1 - cos(theta)) / pow(theta, 2)))) /
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
    return -K(2, 2) * p(1) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2) +
           (-K(1, 1) * p(2) + K(1, 2) * p(1)) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2));
  }
}

inline double dv_dvy(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
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
                         2 * rot(1) * (1 - cos(theta)) / pow(theta, 2)))) /
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
    return -K(1, 2) * p(0) /
               (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2)) +
           K(2, 2) * p(0) *
               (K(1, 1) * T(1) + K(1, 2) * T(2) +
                p(0) * (K(1, 1) * rot(2) - K(1, 2) * rot(1)) +
                p(1) * (K(1, 1) + K(1, 2) * rot(0)) +
                p(2) * (-K(1, 1) * rot(0) + K(1, 2))) /
               pow(K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) +
                       K(2, 2) * p(1) * rot(0) + K(2, 2) * p(2),
                   2);
  }
}

inline double dv_dvz(const Eigen::Vector3d &T, const Eigen::Vector3d &rot,
                     const Eigen::Matrix3d &K, const Eigen::Vector3d &p,
                     const Eigen::Vector3d &x) {
  double theta = std::pow(
      std::pow(rot(0), 2) + std::pow(rot(1), 2) + std::pow(rot(2), 2), 0.5);
  if (std::numeric_limits<double>::epsilon() < rot.dot(rot)) {
    return (p(0) * (K(1, 1) *
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
    return K(1, 1) * p(0) /
           (K(2, 2) * T(2) - K(2, 2) * p(0) * rot(1) + K(2, 2) * p(1) * rot(0) +
            K(2, 2) * p(2));
  }
}

#pragma GCC diagnostic pop

} // namespace optimization