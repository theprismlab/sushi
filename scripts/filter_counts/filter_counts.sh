#!/bin/bash

echo Starting filter_counts...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$SAMPLE_META" ]
then
	echo SAMPLE_META parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$SAMPLE_META" = /* ]]
then
	SAMPLE_META=$(ls $SAMPLE_META)
else
	SAMPLE_META=$BUILD_DIR/$SAMPLE_META
fi

#Enforces abs paths
if [[ "$PRISM_BARCODE_COUNTS" = /* ]]
then
	PRISM_BARCODE_COUNTS=$(ls $PRISM_BARCODE_COUNTS)
else
	PRISM_BARCODE_COUNTS=$BUILD_DIR/$PRISM_BARCODE_COUNTS
fi

echo $CELL_LINE_META

#Enforces abs paths
if [[ "$CELL_LINE_META" = /* ]]
then
	CELL_LINE_META=$(ls $CELL_LINE_META)
else
	CELL_LINE_META=$BUILD_DIR/$CELL_LINE_META
fi

echo $CELL_SET_AND_POOL_META

#Enforces abs paths
if [[ "$CELL_SET_AND_POOL_META" = /* ]]
then
	CELL_SET_AND_POOL_META=$(ls $CELL_SET_AND_POOL_META)
else
	CELL_SET_AND_POOL_META=$BUILD_DIR/$CELL_SET_AND_POOL_META
fi

# Enforces abs paths
if [[ "$SKIPPED_WELLS" = /* ]]
then
  SKIPPED_WELLS=$(ls $SKIPPED_WELLS)
else
  SKIPPED_WELLS=$BUILD_DIR/$SKIPPED_WELLS
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META
echo PRISM_BARCODE_COUNTS is: $PRISM_BARCODE_COUNTS
echo CELL_LINE_META is: $CELL_LINE_META
echo CONTROL_BARCODE_META is: $BUILD_DIR/CB_meta.csv
echo CELL_SET_AND_POOL_META is: $CELL_SET_AND_POOL_META
echo ID_COLS is: $ID_COLS
echo SCREEN is: $SCREEN
echo API_KEY is: $API_KEY
echo FILTER_SKIPPED_WELLS is: $FILTER_SKIPPED_WELLS
echo SKIPPED_WELLS is: $SKIPPED_WELLS
echo REMOVE_DATA is: $REMOVE_DATA

args=(
--prism_barcode_counts "$PRISM_BARCODE_COUNTS"
--sample_meta "$SAMPLE_META"
--cell_line_meta "$CELL_LINE_META"
--CB_meta "$BUILD_DIR/CB_meta.csv"
--cell_set_and_pool_meta "$CELL_SET_AND_POOL_META"
--id_cols "$ID_COLS"
--annot_file "${ANNOTATED_COUNTS}"
--filt_file "${FILTERED_COUNTS}"
--out "$BUILD_DIR"
--rm_data "$REMOVE_DATA"
--screen "$SCREEN"
--api_key "$API_KEY"
--filter_skipped_wells "$FILTER_SKIPPED_WELLS"
--skipped_wells "$SKIPPED_WELLS"
)

echo Rscript filter_counts/filter_counts.R "${args[@]}"

Rscript filter_counts/filter_counts.R "${args[@]}"