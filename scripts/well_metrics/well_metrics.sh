#!/bin/bash

echo Starting well_metrics...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi
echo Build dir is: $BUILD_DIR

#Enforces abs paths

enforce_abs_path() {
  local var_name=$1
  local value=${!var_name}

  if [[ "$value" = /* ]]; then
    eval "$var_name=$(ls "$value")"
  else
    eval "$var_name=$BUILD_DIR/$value"
  fi
  echo "$var_name is: ${!var_name}"
}

enforce_abs_path UNKNOWN_BARCODE_COUNTS
enforce_abs_path PRISM_BARCODE_COUNTS
enforce_abs_path ANNOTATED_COUNTS
enforce_abs_path NORMALIZED_COUNTS
enforce_abs_path CELL_SET_AND_POOL_META
enforce_abs_path CELL_LINE_META
enforce_abs_path SAMPLE_META
enforce_abs_path QC_PARAMS

args=(
--unknown_barcode_counts "${UNKNOWN_BARCODE_COUNTS}"
--prism_barcode_counts "${PRISM_BARCODE_COUNTS}"
--annotated_counts "${ANNOTATED_COUNTS}"
--normalized_counts "${NORMALIZED_COUNTS}"
--cell_set_and_pool_meta "${CELL_SET_AND_POOL_META}"
--cb_meta "${BUILD_DIR}/CB_meta.csv"
--cell_line_meta "${CELL_LINE_META}"
--sample_meta "${SAMPLE_META}"
--qc_params "${QC_PARAMS}"
--id_cols "${ID_COLS}"
--cell_line_cols "${CELL_LINE_COLS}"
--negcon_type "${CTL_TYPES}"
--count_threshold "${COUNT_THRESHOLD}"
--pseudocount "${PSEUDOCOUNT}"
--filter_qc_flags "${FILTER_QC_FLAGS}"
--out "${BUILD_DIR}"
)

echo Rscript well_metrics/well_metrics.R "${args[@]}"
Rscript well_metrics/well_metrics.R "${args[@]}"
