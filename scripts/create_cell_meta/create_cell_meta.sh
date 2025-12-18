#!/bin/bash

echo Starting create_celldb_metadata...

export API_KEY=$(cat /local/jenkins/.clue_api_key)
export API_URL="https://api.clue.io/api/"

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

case "$SEQ_TYPE" in
  HiSeq|NovaSeq)
    INDEX_1="barcode_1"
    INDEX_2="barcode_2"
    BARCODE_SUFFIX="unmapped.1."
    ;;
  MiSeq)
    INDEX_1="_I1_"
    INDEX_2="_I2_"
    BARCODE_SUFFIX="_R1_"
    ;;
  DRAGEN)
    INDEX_1="_I1_"
    INDEX_2="_I2_"
    BARCODE_SUFFIX="_R1_"
    ;;
  *)
	echo "Unknown SEQ_TYPE"
    exit -1
    ;;
esac

if [[ "$SEQ_TYPE" == "NovaSeq" ]]
then
	export REVERSE_INDEX2=TRUE
    echo Sequencing type is NovaSeq
fi

if [[ "$SEQ_TYPE" == "DRAGEN" ]]
then
	export REVERSE_INDEX2=TRUE
    echo Sequencing type is DRAGEN
fi

if [ -z "$INDEX_1" ]
then
	echo INDEX_1 parameter empty
    exit -1
fi

if [ -z "$INDEX_2" ]
then
	echo INDEX_2 parameter empty
    exit -1
fi

if [ -z "$BARCODE_SUFFIX" ]
then
	echo BARCODE_SUFFIX parameter empty
    exit -1
fi


#Enforces abs paths
if [[ "$SAMPLE_META" = /* ]]
then
	SAMPLE_META=$(ls $SAMPLE_META)
else
	SAMPLE_META=$BUILD_DIR/$SAMPLE_META
fi

echo Build dir is: $BUILD_DIR

PROJECT_DIR=$(dirname "$BUILD_DIR")
PROJECT_CODE=$(basename "$PROJECT_DIR")

echo Project Code: $PROJECT_CODE

# Define an array of parameters and their corresponding values
parameters=(
  "BUILD_DIR:$BUILD_DIR"
  "SAMPLE_META:$SAMPLE_META"
  "SEQ_TYPE:$SEQ_TYPE"
  "INDEX_1:$INDEX_1"
  "INDEX_2:$INDEX_2"
  "BARCODE_SUFFIX:$BARCODE_SUFFIX"
  "RUN_NORM:$RUN_NORM"
  "ID_COLS:$ID_COLS"
  "SIG_COLS:$SIG_COLS"
  "SEQUENCING_INDEX_COLS:$SEQUENCING_INDEX_COLS"
  "CONTROL_COLS:$CONTROL_COLS"
  "COUNT_THRESHOLD:$COUNT_THRESHOLD"
  "PSEUDOCOUNT:$PSEUDOCOUNT"
  "COUNT_COL_NAME:$COUNT_COL_NAME"
  "CONTROL_BARCODE_META:$CONTROL_BARCODE_META"
  "BUILD_NAME:$BUILD_NAME"
  "API_KEY:$API_KEY"
  "CELL_LINE_COLS:$CELL_LINE_COLS"
)

echo $R_LIBS

args=(
--sample_meta "$SAMPLE_META"
--out "$BUILD_DIR"
--cb_ladder "$CONTROL_BARCODE_META"
--api_key "$API_KEY"
--cell_line_cols "$CELL_LINE_COLS"
)

echo Rscript create_cell_meta/create_cell_meta.R "${args[@]}"

Rscript create_cell_meta/create_cell_meta.R "${args[@]}"
