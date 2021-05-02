#!/bin/bash

SCRIPT=$(readlink -f "$0")
CUR_DIR=$(dirname "$SCRIPT")

source ${CUR_DIR}/../env.sh

${CUR_DIR}/../install/bin/appl_bundle_adjustment_gradient_descent \
    --camera_extrinsic_file_path=/home/koichi/private_data/sample_data/camera_extrinsics.txt \
    --camera_intrinsic_file_path=/home/koichi/private_data/sample_data/camera_intrinsics.txt \
    --image_location_file_path=/home/koichi/private_data/sample_data/image_points2d.txt \
    --points_location_path=/home/koichi/private_data/sample_data/points3d.txt \
    --result_save_dir=/home/koichi/private_data/sample_data
