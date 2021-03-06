#!/bin/bash

source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/build_env.sh"

# Build setting.
CMAKE_BUILD_TYPE=Release
#CMAKE_BUILD_TYPE=Debug
COMPILER_WARNING="-Wall -Wextra"

BUILD_TEST_PKG=TRUE

if [ ! -e ./build ]; then
  mkdir build
fi

cd build

${CMAKE_BIN} \
  -D CMAKE_C_COMPILER=gcc-9 \
  -D CMAKE_CXX_COMPILER=g++-9 \
  -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -D CMAKE_VERBOSE_MAKEFILE=FALSE \
  -D COMPILER_WARNING=${COMPILER_WARNING} \
  -D BUILD_TEST_PKG=${BUILD_TEST_PKG} \
  -D gflags_DIR=${gflags_DIR} \
  -D Glog_DIR=${glog_DIR} \
  -D GTEST_ROOT=${gtest_DIR} \
  -D OpenCV_DIR=${OpenCV_DIR} \
  -D Ceres_DIR=${Ceres_DIR} \
  ../

make install -j4

cd ${PROJECT_DIR}