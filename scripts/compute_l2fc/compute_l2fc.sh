#!/bin/bash

echo Starting compute_l2fc...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$NORMALIZED_COUNTS" ]
then
	echo NORMALIZED_COUNTS parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$NORMALIZED_COUNTS" = /* ]]
then
	NORMALIZED_COUNTS=$(ls $NORMALIZED_COUNTS)
else
	NORMALIZED_COUNTS=$BUILD_DIR/$NORMALIZED_COUNTS
fi

# Construct path to plate+cell qc table
PLATE_CELL_QC_TABLE="$BUILD_DIR"/qc_tables/plate_cell_qc_table.csv

echo Build dir is: $BUILD_DIR
echo NORMALIZED_COUNTS is: $NORMALIZED_COUNTS

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--output_file "${LFC}"
--out "$BUILD_DIR"
--control_type "$CTL_TYPES"
--count_col_name "$COUNT_COL_NAME"
--sig_cols "$SIG_COLS"
--bio_rep_col "${BIO_REP_COL}"
--ctrl_cols "$CONTROL_COLS"
--cell_line_cols "$CELL_LINE_COLS"
--qc_path "$PLATE_CELL_QC_TABLE"
--filter_failed_lines "$FILTER_FAILED_LINES"
)

echo Rscript compute_l2fc/compute_l2fc.R "${args[@]}"
Rscript compute_l2fc/compute_l2fc.R "${args[@]}"
