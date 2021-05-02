
cmake_minimum_required(VERSION 3.8)
project(sfm_exp)

cmake_policy(SET CMP0074 NEW)
set(Boost_USE_RELEASE_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

# CMake Compile Option
add_compile_options(${COMPILE_OPTION})

# CXX Version.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_WARNING}")

# Care for boost.
add_definitions(-DBOOST_ALL_NO_LIB)
find_package(Boost REQUIRED COMPONENTS system thread)

# External library management.
find_package(Eigen3 REQUIRED)
find_package(Gflags REQUIRED)
find_package(Glog REQUIRED)
find_package(GTest REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(PCL REQUIRED)

# Module directory.
add_subdirectory(src/utility)
add_subdirectory(src/core)

# Application directory. All executables contained in this directory.
add_subdirectory(src/applications)