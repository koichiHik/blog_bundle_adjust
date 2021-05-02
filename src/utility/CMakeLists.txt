
cmake_minimum_required(VERSION 3.8)
project(sfm_exp)

set(LIBRARY_NAME utility)
################################
# Description for library.
################################
set(SRC_FILES
    src/file_utility.cc
    src/visualize.cc
)
set(PUBLIC_HEADERS
    "include/file_utility.h"
    "include/visualize.h"
)

add_library(${LIBRARY_NAME} SHARED ${SRC_FILES})
set_target_properties(${LIBRARY_NAME} PROPERTIES PUBLIC_HEADER 
    "${PUBLIC_HEADERS}"
)
target_include_directories(${LIBRARY_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
    include
)
target_include_directories(${LIBRARY_NAME} PRIVATE
    ${PCL_INCLUDE_DIRS}
)
target_link_libraries(${LIBRARY_NAME} PUBLIC
    optimization
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
    ${PCL_VISUALIZATION_LIBRARIES}
)

install(TARGETS ${LIBRARY_NAME} 
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

################################
# Description for testing.
################################
# Set whole core directory as include directory.
#include_directories(${GTEST_INCLUDE_DIRS})

# Test for functions in unnamed namespacd in depth_image_util.cc
#add_executable(core_unit_test
#    src/geometry.cc
#    test/test_rq_decomposition.cc
#)
#target_include_directories(core_unit_test PUBLIC 
#    ${EIGEN3_INCLUDE_DIR}
#    include
#    src
#)
#target_link_libraries(core_unit_test PRIVATE
#    GTest::GTest 
#    GTest::Main
#    glog::glog
#)
#target_compile_definitions(core_unit_test PRIVATE CANLIB_EXPORT=1)
#
# Google Testの各テストケースごとにCTestのテストを作成する
#gtest_add_tests(TARGET test_depth_image_util_unnamed)

# CTestのテストをひとつだけ作成する
#add_test(NAME core_unit_test 
#    COMMAND core_unit_test)