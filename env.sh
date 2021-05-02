#!/bin/bash

SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

DEP_DIR="/home/koichi/private_ws/3rdParty/install"
gflags_LIB_DIR="${DEP_DIR}/gflags/lib"
glog_LIB_DIR="${DEP_DIR}/glog/lib"
gtest_LIB_DIR="${DEP_DIR}/gtest/lib"
OpenCV_LIB_DIR="${DEP_DIR}/opencv/lib"

LD_LIBRARY_PATH=${gflags_LIB_DIR}:${glog_LIB_DIR}:${gtest_LIB_DIR}:${OpenCV_LIB_DIR}:$LD_LIBRARY_PATH

INTERNAL_LIB_DIR="${SCRIPT_DIR}/../install/lib"

LD_LIBRARY_PATH=${INTERNAL_LIB_DIR}:$LD_LIBRARY_PATH
