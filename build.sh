#!/bin/bash

# Install directory is relative to this script.
TOP_DIR=$(cd $(dirname $0); pwd)
INSTALL_DIR="${TOP_DIR}/install"
DEP_DIR="/home/koichi/private_ws/3rdParty/install"

# CMake path.
CMAKE_BIN="${DEP_DIR}/cmake/bin/cmake"

# Third party library.
gflags_DIR="${DEP_DIR}/gflags/lib/cmake/gflags"
glog_DIR="${DEP_DIR}/glog/lib/cmake/glog"
gtest_DIR="${DEP_DIR}/gtest/lib"
OpenCV_DIR="${DEP_DIR}/opencv/lib/cmake/opencv4"
cereal_DIR="${DEP_DIR}/cereal/share/cmake/cereal"
#PCL_DIR="${DEP_DIR}/pcl/share/pcl-1.11"
leveldb_DIR="${DEP_DIR}/leveldb/lib/cmake/leveldb"
Ceres_DIR="/usr/lib/cmake/Ceres"

# Third party library. (Non CMAKE)
vlfeat_INCLUDE_DIR="${DEP_DIR}/vlfeat/include"
vlfeat_LIBRARIES="${DEP_DIR}/vlfeat/lib/libvl.so"

# Build setting.
CMAKE_BUILD_TYPE=Release
#CMAKE_BUILD_TYPE=Debug
COMPILER_WARNING="-Wall -Wextra"

BUILD_TEST_PKG=TRUE

if [ ! -e ./build ]; then
  mkdir build
fi

cd build

#  -D PCL_DIR=${PCL_DIR} \


${CMAKE_BIN} \
  -D CMAKE_C_COMPILER=gcc-9 \
  -D CMAKE_CXX_COMPILER=g++-9 \
  -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -D CMAKE_CONFIGURATION_TYPES=${CMAKE_BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -D CMAKE_VERBOSE_MAKEFILE=FALSE \
  -D COMPILER_WARNING=${COMPILER_WARNING} \
  -D BUILD_TEST_PKG=${BUILD_TEST_PKG} \
  -D gflags_DIR=${gflags_DIR} \
  -D Glog_DIR=${glog_DIR} \
  -D GTEST_ROOT=${gtest_DIR} \
  -D OpenCV_DIR=${OpenCV_DIR} \
  -D cereal_DIR=${cereal_DIR} \
  -D vlfeat_INCLUDE_DIR=${vlfeat_INCLUDE_DIR} \
  -D vlfeat_LIBRARIES=${vlfeat_LIBRARIES} \
  -D leveldb_DIR=${leveldb_DIR} \
  -D Ceres_DIR=${Ceres_DIR} \
  ../

make install -j4 && bash ${TOP_DIR}/script/run_unit_test.sh

cd ${TOP_DIR}