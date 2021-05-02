#!/bin/bash

SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

DEP_DIR="/home/koichi/private_ws/3rdParty/install"
gflags_LIB_DIR="${DEP_DIR}/gflags/lib"
glog_LIB_DIR="${DEP_DIR}/glog/lib"
gtest_LIB_DIR="${DEP_DIR}/gtest/lib"
OpenCV_LIB_DIR="${DEP_DIR}/opencv/lib"
vl_LIB_DIR="${DEP_DIR}/vlfeat/lib"
PCL_LIB_DIR="${DEP_DIR}/pcl/lib"
leveldb_LIB_DIR="${DEP_DIR}/leveldb/lib"

LD_LIBRARY_PATH=${gflags_LIB_DIR}:${glog_LIB_DIR}:${gtest_LIB_DIR}:${OpenCV_LIB_DIR}:${vl_LIB_DIR}:${PCL_LIB_DIR}:${leveldb_LIB_DIR}:$LD_LIBRARY_PATH

INTERNAL_LIB_DIR="${SCRIPT_DIR}/../install/lib"

LD_LIBRARY_PATH=${INTERNAL_LIB_DIR}:$LD_LIBRARY_PATH