
cmake_minimum_required(VERSION 3.8)
project(sfm_exp)

#############################
# Construct Artificial Data.
#############################
set(APP_NAME appl_create_artificial_data)
add_executable(${APP_NAME}
    appl_create_artificial_data.cc)
target_include_directories(${APP_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
)
target_link_libraries(${APP_NAME} PUBLIC
    optimization
    utility
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
    ${PCL_VISUALIZATION_LIBRARIES}
)
install(TARGETS ${APP_NAME} 
    RUNTIME DESTINATION bin)

#############################
# Application with gradient descent.
#############################
set(APP_NAME appl_bundle_adjustment_gradient_descent)
add_executable(${APP_NAME}
    appl_bundle_adjustment_gradient_descent.cc)
target_include_directories(${APP_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
)
target_link_libraries(${APP_NAME} PUBLIC
    optimization
    utility
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
    ${PCL_VISUALIZATION_LIBRARIES}
)
install(TARGETS ${APP_NAME} 
    RUNTIME DESTINATION bin)

#############################
# Application with pure newton.
#############################
set(APP_NAME appl_bundle_adjustment_pure_newton)
add_executable(${APP_NAME}
    appl_bundle_adjustment_pure_newton.cc)
target_include_directories(${APP_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
)
target_link_libraries(${APP_NAME} PUBLIC
    optimization
    utility
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
    ${PCL_VISUALIZATION_LIBRARIES}
)
install(TARGETS ${APP_NAME} 
    RUNTIME DESTINATION bin)

#############################
# Application with ceres solver.
#############################
set(APP_NAME appl_bundle_adjustment_ceres_solver)
add_executable(${APP_NAME}
    appl_bundle_adjustment_ceres_solver.cc)
target_include_directories(${APP_NAME} PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
)
target_link_libraries(${APP_NAME} PUBLIC
    optimization
    utility
    glog::glog
    ${GFLAGS_LIBRARIES}
    ${OpenCV_LIBS}
    ${PCL_VISUALIZATION_LIBRARIES}
)
install(TARGETS ${APP_NAME} 
    RUNTIME DESTINATION bin)