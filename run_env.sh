#!/bin/bash

# Directory of this script.
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

THIRD_PARTY_DIR="${PROJECT_DIR}/../3rdparty/install"

# Library
gflags_LIB_DIR="${THIRD_PARTY_DIR}/gflags/lib"
glog_LIB_DIR="${THIRD_PARTY_DIR}/glog/lib"
gtest_LIB_DIR="${THIRD_PARTY_DIR}/googletest/lib"
OpenCV_LIB_DIR="${THIRD_PARTY_DIR}/opencv/lib"
Ceres_LIB_DIR="/usr/lib"

LD_LIBRARY_PATH=${gflags_LIB_DIR}:${glog_LIB_DIR}:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=${gtest_LIB_DIR}:${OpenCV_LIB_DIR}:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=${Ceres_LIB_DIR}:$LD_LIBRARY_PATH

INTERNAL_LIB_DIR="${PROJECT_DIR}/install/lib"
LD_LIBRARY_PATH=${INTERNAL_LIB_DIR}:$LD_LIBRARY_PATH
