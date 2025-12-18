library(argparse)
library(tidyverse)

source("utils/kitchen_utensils.R")
source("bias_correction/bias_correction_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values.")
parser$add_argument("--sig_cols", default = "pert_name,pert_dose,pert_dose_unit,day",
                    help = "Columns used to generate signature ids")
parser$add_argument("--bio_rep_col", default = "", help = "Column that identifies the biological replicate.")
parser$add_argument("--l2fc_col", default = "l2fc", help = "Column containing log2 fold change values to be corrected.")
parser$add_argument("--growth_pattern_col", default = "growth_condition",
                    help = "Column containing growth pattern annotations.")
args = parser$parse_args()

# Read in files and set up any args
l2fc = read_data_table(args$l2fc)
sig_cols = unlist(strsplit(args$sig_cols, ","))
bio_rep_col = args$bio_rep_col
l2fc_col = args$l2fc_col
growth_pattern_col = args$growth_pattern_col

# Combine sig_cols and bio_rep
if (bio_rep_col == "") {
  message("No bio_rep column specified. Assuming that there are NO biological replicates.")
  sig_bio_rep_cols = sig_cols
} else {
  sig_bio_rep_cols = unique(c(sig_cols, bio_rep_col))
}

# Make sure cell_set is not in sig_bio_rep_cols
if ("cell_set" %in% sig_bio_rep_cols) {
  message("Detecting cell_set in sig_bio_cols Dropping this column from the list.")
  sig_bio_rep_cols = sig_bio_rep_cols[sig_bio_rep_cols != "cell_set"]
}

# Throw an error if a cell line does not have a growth annotation
if (any(unique(l2fc[[growth_pattern_col]]) %in% c(NA, "", " ", "NA"))) {
  stop("Detecting empty or missing values in the growth pattern column - ", growth_pattern_col)
}

# Correct l2fcs by regressing out cell line growth patterns
corrected_l2fc = l2fc |>
  dplyr::mutate(negcon_log2_norm_n = log2(control_median_normalized_n)) |>
  dplyr::group_split(dplyr::across(tidyselect::all_of(sig_bio_rep_cols))) |>
  lapply(apply_bias_correction,
         raw_l2fc_col = l2fc_col,
         growth_pattern_col = growth_pattern_col,
         cell_set_col = "cell_set") |>
  dplyr::bind_rows() |>
  dplyr::ungroup()

# Check the column has been generated
if ("l2fc_uncorrected" %in% colnames(corrected_l2fc)) {
  if (!NA %in% unique(corrected_l2fc$l2fc_uncorrected)) {
    message("Success: l2fc_uncorrected column was generated.")
  } else {
    message("Warning: l2fc_uncorrected column was generated but with NA values.")
  }
  # Write out
  message("Writing corrected l2fc value to ", args$l2fc)
  write_out_table(corrected_l2fc, args$l2fc)
  
} else {
  stop("Error: l2fc_uncorrected was not generated.")
}

# Write out
message("Writing corrected l2fc value to ", args$l2fc)
write_out_table(corrected_l2fc, args$l2fc)

# Check file creation
check_file_exists(args$l2fc)