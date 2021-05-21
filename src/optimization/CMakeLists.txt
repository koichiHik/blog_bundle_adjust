
cmake_minimum_required(VERSION 3.8)
project(sfm_exp)

set(LIBRARY_NAME optimization)
################################
# Description for library.
################################
set(SRC_FILES
    src/geometry.cc
    #src/bundle_adjust.cc
    #src/matrix_computation.cc
    src/error_metric.cc
    src/gradient_descent.cc
    src/pure_newton.cc
    src/differential_operators.cc
    src/compute_hessian.cc
    src/line_searches.cc
    src/ceres_solver.cc
)
set(PUBLIC_HEADERS
    "include/geometry.h"
    "include/common_def.h"
)

add_library(${LIBRARY_NAME} SHARED ${SRC_FILES})
set_target_properties(${LIBRARY_NAME} PROPERTIES PUBLIC_HEADER 
    "${PUBLIC_HEADERS}"
)
target_include_directories(${LIBRARY_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    include
)
target_link_libraries(${LIBRARY_NAME} PUBLIC
    ceres
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
)

install(TARGETS ${LIBRARY_NAME} 
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

################################
# Description for testing.
################################
# Set whole core directory as include directory.
include_directories(${GTEST_INCLUDE_DIRS})

# Test for functions in unnamed namespacd in depth_image_util.cc
add_executable(optimization_unit_test
    src/geometry.cc
    test/test_rq_decomposition.cc
)
target_include_directories(optimization_unit_test PUBLIC 
    ${EIGEN3_INCLUDE_DIR}
    include
    src
)
target_link_libraries(optimization_unit_test PRIVATE
    GTest::GTest 
    GTest::Main
    glog::glog
    ${OpenCV_LIBS}
)
target_compile_definitions(optimization_unit_test PRIVATE CANLIB_EXPORT=1)

# Google Testの各テストケースごとにCTestのテストを作成する
#gtest_add_tests(TARGET test_depth_image_util_unnamed)

# CTestのテストをひとつだけ作成する
add_test(NAME optimization_unit_test 
    COMMAND optimization_unit_test)