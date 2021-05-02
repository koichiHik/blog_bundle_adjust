#!/bin/bash

SCRIPT=$(readlink -f "$0")
CUR_DIR=$(dirname "$SCRIPT")

source ${CUR_DIR}/../env.sh

${CUR_DIR}/../install/bin/appl_compute_track_locations \
    --camera_file_path=/home/koichi/private_data/oxford_vis/dinosaur/cameras.txt \
    --image_location_file_path=/home/koichi/private_data/oxford_vis/dinosaur/image_points2d.txt \
    --points_location_path=/home/koichi/private_data/oxford_vis/dinosaur/points3d.txt