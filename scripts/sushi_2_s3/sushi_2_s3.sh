#!/bin/bash

echo Starting sushi_2_s3...

if [ -z "$BUILD_DIR" ]; then
    echo "BUILD_DIR not specified"
    exit -1
fi

echo "Build dir is: $BUILD_DIR"

S3_BUCKET="macchiato.clue.io"
if [ "$TEST_DATASET" = "true" ]; then
    S3_BUCKET="test-macchiato.clue.io"
    echo "TEST_DATASET enabled — uploading to test bucket: $S3_BUCKET"
fi

args=(
--build_path "$BUILD_DIR"
--s3_bucket "$S3_BUCKET"
)

echo python3 sushi_2_s3/sushi_2_s3.py "${args[@]}"

python3 sushi_2_s3/sushi_2_s3.py "${args[@]}"
