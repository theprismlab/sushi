library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
library(rhdf5)

source("utils/kitchen_utensils.R")
source("compute_l2fc/compute_l2fc_functions.R")
source("synergy/synergy_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--wkdir", default = getwd(), help = "Working directory")
parser$add_argument("--normalized_counts", default = "normalized_counts.csv",
                    help = "Path to file containing normalized counts")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values")
# Used to create mock DMSO l2fcs and to restructure l2fcs
parser$add_argument("--cell_line_cols", default = "pool_id,depmap_id,lua",
                    help = "Column names that describe a cell line")
# Used to create mock DMSO l2fcs
parser$add_argument("--ctrl_cols", default = "cell_set,day",
                    help = "Column names that describe the control conditions")
# Used to create mock DMSO l2fcs and to restructure l2fcs
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "Column names that describe a unique profile")
# Used to restructure l2fcs
parser$add_argument("--combination_col", default = "is_combination", help = "Column name used to identify combinations")
# Used to create mock DMSO l2fcs
parser$add_argument("--count_col_name", default = "log2_normalized_n",
                    help = "Column name containing counts to be used to calculate l2fcs")
# Used to restructure l2fcs and to calculate synergy
parser$add_argument("--synergy_l2fc_col", default = "median_l2fc_uncorrected",
                    help = "Column name containing l2fc values")
# Used to create mock DMSO l2fcs
parser$add_argument("--n_samples", default = 10000, help = "Size of the resampling")
# Used to create mock DMSO l2fcs
parser$add_argument("--negcon_type", default = "ctl_vehicle",
                    help = "String in pert_type that identifies the negative control")
# Name of synergy output
parser$add_argument("--output_file", default = "synergy_scores.csv", help = "File to output synergy scores")
parser$add_argument("-o", "--out", default = "", help = "Output path, defaults to working directory")

args = parser$parse_args()

# Read in files and set up parameters ----
normalized_counts = read_data_table(args$normalized_counts)
cps_l2fc = read_data_table(args$l2fc)

args$cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
args$ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
args$sig_cols = unlist(strsplit(args$sig_cols, ","))
args$n_samples = as.numeric(args$n_samples)

# Check if combination_col is present ----
if (!args$combination_col %in% colnames(cps_l2fc)) {
  stop("The combination column is not detected!")
}

# Create DMSO L2FC for synergy pvalue calculations ----
# DMSO L2FCs will be resampled to generate a null distribution for each cell line + pcr_plate

# Pull out negcon wells from norm counts and duplicate with different pert_type
dmso_norm_ctrl = normalized_counts[pert_type == args$negcon_type, ]
dmso_norm_trt = dmso_norm_ctrl
dmso_norm_trt$pert_type = "trt_cp"
dmso_norm_input = data.table::rbindlist(list(dmso_norm_trt, dmso_norm_ctrl), use.names = TRUE)

# Median center dmso values for each pcr plate
dmso_norm_input[, log2_normalized_n := log2_normalized_n - median(log2_normalized_n), by = pcr_plate]

# Call compute_l2fc to get l2fc values
dmso_l2fc = compute_l2fc(normalized_counts = dmso_norm_input,
                         cell_line_cols = args$cell_line_cols,
                         control_type = args$negcon_type,
                         count_col_name = args$count_col_name,
                         ctrl_cols = args$ctrl_cols,
                         sig_cols = args$sig_cols)

# Resample DMSO L2FC to create mock synergy values ----
group_name_cols = unique(c(args$cell_line_cols, args$ctrl_cols)) # Also used later in pval calc
# HOTFIX: Remove pcr_plate and replicate from group_name_cols if present (not in l2fc file)
group_name_cols = group_name_cols[!group_name_cols %in% c("pcr_plate", "replicate_plate")]

create_dmso_synergy_hdf5(dmso_l2fc = dmso_l2fc,
                         group_name_cols = group_name_cols,
                         path = file.path(args$out, "dmso_synergy.h5"),
                         n_samples =  args$n_samples)

# Restructure L2FCs with single agent L2FCs and calculate synergy ----
# Create vector of new column names
names_prefix = c("pert", "pert2", "combo")
new_names = paste(names_prefix, args$synergy_l2fc_col, sep = "_")

# Pull out columns that describe each perturbation
ignore_cols = c("pert_type", "pert_plate", "pert_vehicle") # common columns with "pert" to be ignored
pert1_cols = sort(args$sig_cols[grepl(paste0(names_prefix[1], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols])
pert2_cols = sort(args$sig_cols[grepl(paste0(names_prefix[2], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols])

# Pull out other columns to join on - these are sig_cols that don't describe the perturbations
filt_sig_cols = args$sig_cols[!args$sig_cols %in% c(pert1_cols, pert2_cols, args$combination_col)]

# Call restructure to "pivot" single agent l2fcs to additional columns for combination data
trt_synergy = restructure_l2fc(cps_l2fc = cps_l2fc,
                               join_cols = unique(c(args$cell_line_cols, filt_sig_cols)),
                               pert_cols_list = list(pert1_cols, pert2_cols),
                               l2fc_col = args$synergy_l2fc_col,
                               new_l2fc_col_names = new_names,
                               combination_col = args$combination_col)

# Add synergy scores by comparing the single agent l2fc columns to the combination l2fc column
calculate_synergy(restructured_l2fc = trt_synergy, l2fc_cols = new_names)
print("Calculated synergies.")

# Pull out pvalues ----
# Open hdf5 connection
dmso_synergy = rhdf5::H5Fopen(file.path(args$out, "dmso_synergy.h5"))
existing_groups = rhdf5::h5ls(file.path(args$out, "dmso_synergy.h5"))

# Create group name and filter for only cls that were sampled
trt_synergy[, group_name := do.call(paste, c(.SD, sep = "__")), .SDcols = group_name_cols]
# Filter out names not in hf file
trt_synergy = trt_synergy[group_name %in% existing_groups$name, ]

# Calculate emperical pvalues using mapply
trt_synergy[, p_val_emp := mapply(get_pvalue, group_name, synergy, MoreArgs = list(h5_file = dmso_synergy))]
rhdf5::h5closeAll()
print("Added empirical pvalues to synergies.")

# Calculate qvalues
trt_synergy[, q_val := p.adjust(p_val_emp, method = "BH"), by = c(args$sig_cols)]

# Count number of significant cell lines for a combination profile
sig_cols_wo_dose = args$sig_cols[!grepl("_dose", args$sig_cols)]
trt_synergy[, num_sig_profiles := sum(q_val < 0.005), by = c(args$cell_line_cols, sig_cols_wo_dose)]

# Write out file ----
outpath = file.path(args$out, args$output_file)
print(paste0("Writing out synergy file to ", outpath))
write_out_table(table = trt_synergy, path = outpath)

# End ----