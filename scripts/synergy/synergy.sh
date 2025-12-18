
echo Starting synergy module...

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
enforce_abs_path COLLAPSED_LFC

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--l2fc "$COLLAPSED_LFC"
--cell_line_cols "$CELL_LINE_COLS"
--ctrl_cols "$CONTROL_COLS"
--sig_cols "$SIG_COLS"
--combination_col "$COMBINATION_COL"
--count_col_name "$COUNT_COL_NAME" 
--synergy_l2fc_col "$SYNERGY_L2FC_COL"
--n_samples "$N_SAMPLES"
--negcon_type "$CTL_TYPES"
--output_file "${SYNERGY_FILE}"
--out "$BUILD_DIR"
)

echo Rscript synergy/synergy.R "${args[@]}"
Rscript synergy/synergy.R "${args[@]}"
