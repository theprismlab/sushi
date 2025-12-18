echo Starting bias correction module ....

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi
echo Build dir is: $BUILD_DIR

# Function to help build paths
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

enforce_abs_path LFC

args=(
--l2fc "${LFC}"
--sig_cols "${SIG_COLS}"
--bio_rep_col "${BIO_REP_COL}"
--l2fc_col "${L2FC_COLUMN}"
--growth_pattern_col "${GROWTH_PATTERN_COL}"
)

echo Rscript bias_correction/bias_correction.R "${args[@]}"
Rscript bias_correction/bias_correction.R "${args[@]}"