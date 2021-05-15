#!/bin/bash

# Directory of this script.
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

THIRD_PARTY_DIR="${PROJECT_DIR}/../3rdparty/install"

# Install Directory
INSTALL_DIR="${PROJECT_DIR}/install"

# CMake
CMAKE_BIN="${THIRD_PARTY_DIR}/cmake/bin/cmake"

# Library
gflags_DIR="${THIRD_PARTY_DIR}/gflags/lib/cmake/gflags"
glog_DIR="${THIRD_PARTY_DIR}/glog/lib/cmake/glog"
gtest_DIR="${THIRD_PARTY_DIR}/googletest"
OpenCV_DIR="${THIRD_PARTY_DIR}/opencv/lib/cmake/opencv4"
Ceres_DIR="/usr/lib/cmake/Ceres"