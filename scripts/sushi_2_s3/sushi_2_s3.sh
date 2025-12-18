#!/bin/bash

echo Starting sushi_2_s3...

if [ -z "$BUILD_DIR" ]; then
    echo "BUILD_DIR not specified"
    exit -1
fi

echo "Build dir is: $BUILD_DIR"

args=(
--build_path "$BUILD_DIR"
--s3_bucket "macchiato.clue.io"
)

echo python3 sushi_2_s3/sushi_2_s3.py "${args[@]}"

python3 sushi_2_s3/sushi_2_s3.py "${args[@]}"
