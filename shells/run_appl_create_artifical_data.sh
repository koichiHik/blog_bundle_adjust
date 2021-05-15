#!/bin/bash

source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/../run_env.sh"

${PROJECT_DIR}/install/bin/appl_create_artificial_data \
    --camera_extrinsic_file_path=/home/koichi/private_data/sample_data/camera_extrinsics.txt \
    --camera_intrinsic_file_path=/home/koichi/private_data/sample_data/camera_intrinsics.txt \
    --image_location_file_path=/home/koichi/private_data/sample_data/image_points2d.txt \
    --points_location_path=/home/koichi/private_data/sample_data/points3d.txt