#!/bin/bash

echo Starting generate_qc_tables...

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

enforce_abs_path NORMALIZED_COUNTS
enforce_abs_path QC_PARAMS

args=(
--normalized_counts "${NORMALIZED_COUNTS}"
--qc_params "${QC_PARAMS}"
--id_cols "${ID_COLS}"
--cell_line_cols "${CELL_LINE_COLS}"
--sig_cols "${SIG_COLS}"
--pseudocount "${PSEUDOCOUNT}"
--negcon_type "${CTL_TYPES}"
--poscon_type "${POSCON_TYPE}"
--out "${BUILD_DIR}"
)

echo Rscript qc_tables/qc_tables.R "${args[@]}"
Rscript qc_tables/qc_tables.R "${args[@]}"
