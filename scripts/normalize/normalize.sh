#!/bin/bash

echo Starting CBnormalize...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$FILTERED_COUNTS" ]
then
	echo FILTERED_COUNTS parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$FILTERED_COUNTS" = /* ]]
then
	FILTERED_COUNTS=$(ls $FILTERED_COUNTS)
else
	FILTERED_COUNTS=$BUILD_DIR/$FILTERED_COUNTS
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META


echo $RUN_NORM


echo "Running normalization module"

echo Rscript normalize/normalize.R -c $FILTERED_COUNTS	\
--CB_meta $BUILD_DIR/CB_meta.csv \
--norm_method ${NORM_METHOD} \
--pseudocount $PSEUDOCOUNT \
--id_cols $ID_COLS \
--read_detection_limit $READ_DETECTION_LIMIT \
--negcon_cols $CONTROL_COLS \
--negcon_type $CTL_TYPES \
--output_file $NORMALIZED_COUNTS \
--out $BUILD_DIR

Rscript normalize/normalize.R -c $FILTERED_COUNTS	\
--CB_meta $BUILD_DIR/CB_meta.csv \
--norm_method ${NORM_METHOD} \
--pseudocount $PSEUDOCOUNT \
--id_cols $ID_COLS \
--read_detection_limit $READ_DETECTION_LIMIT \
--negcon_cols $CONTROL_COLS \
--negcon_type $CTL_TYPES \
--output_file $NORMALIZED_COUNTS \
--out $BUILD_DIR

