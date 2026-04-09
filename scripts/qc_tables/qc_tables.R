options(cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(dplyr)
})
source("qc_tables/qc_tables_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser = ArgumentParser()
parser$add_argument("--normalized_counts", default = "normalized_counts.csv", help = "normalized counts file")
parser$add_argument("--qc_params", default = "qc_params.json", help = "File containing QC parameters")
parser$add_argument("--cell_line_cols", default = "pool_id,depmap_id,lua")
parser$add_argument("--sig_cols", default = "")
parser$add_argument("--pcr_plate_col", default = "pcr_plate")
parser$add_argument("--pert_plate_col", default = "pert_plate")
parser$add_argument("--pseudocount", default = 20)
parser$add_argument("-n", "--negcon_type", default = "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default = "trt_poscon")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")
args = parser$parse_args()

# Read in files and set up parameters ----
message("Reading in normalized_counts from", args$normalized_counts, "...")
normalized_counts = read_data_table(args$normalized_counts)
message("Reading in qc_thresholds from ", args$qc_params, "...")
thresholds = load_thresholds_from_json(args$qc_params)

# Create a qc output directory if one does not exist
if (!dir.exists(file.path(args$out, "qc_tables"))) {
  dir.create(file.path(args$out, "qc_tables"))
}

# Set up parameters
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
sig_cols = unlist(strsplit(args$sig_cols, ","))
pcr_plate_col = args$pcr_plate_col
pert_plate_col = args$pert_plate_col
pseudocount = as.numeric(args$pseudocount)
negcon = args$negcon_type
poscon = args$poscon_type
contains_poscon = any(normalized_counts$pert_type == args$poscon_type)

# Filter out control barcodes from normalized counts
normalized_counts_rm_cbc = filter_control_barcodes(normalized_counts)

# Outlier pools ----
# Identify outlier pools
outlier_pools = get_outlier_pools(normalized_counts_rm_cbc,
                                  negcon = negcon,
                                  id_cols = c("pcr_plate", "pcr_well"),
                                  pert_plate_col = pert_plate_col,
                                  pool_cols = c("cell_set", "pool_id"),
                                  cell_line_cols = cell_line_cols)

# Filter out poor pools
flagged_pools = outlier_pools |> dplyr::filter(!is.na(well_flag))

norm_counts_filt_pools = dplyr::anti_join(normalized_counts_rm_cbc, flagged_pools,
                                          by = c(id_cols, "pool_id"))

# Write out files
pool_well_qc_table_outpath = file.path(args$out, "qc_tables", "pool_well_qc_table.csv")
message("Writing out pool_well_qc_table to ", pool_well_qc_table_outpath)
write_out_table(table = outlier_pools, path = pool_well_qc_table_outpath)
check_file_exists(pool_well_qc_table_outpath)

# Write out flagged pools
pool_well_qc_flags_outpath = file.path(args$out, "qc_tables", "pool_well_qc_flags.csv")
message("Writing out pool_well_qc_flags to ", pool_well_qc_flags_outpath)
write_out_table(table = flagged_pools, path = pool_well_qc_flags_outpath)
check_file_exists(pool_well_qc_flags_outpath)

# Plate cell line QCs ----
# Generate plate cell QCs
plate_cell_table = generate_plate_cell_table(
  normalized_counts = normalized_counts_rm_cbc,
  cell_line_cols = cell_line_cols,
  sig_cols = sig_cols,
  pcr_plate_col = pcr_plate_col,
  pert_plate_col = pert_plate_col,
  optional_cols = c("project_code", "day"),
  pseudocount = pseudocount,
  contains_poscon = contains_poscon,
  poscon = poscon, negcon = negcon,
  nc_variability_threshold = thresholds$nc_variability_threshold,
  error_rate_threshold = thresholds$error_rate_threshold,
  pc_viability_threshold = thresholds$pc_viability_threshold,
  nc_raw_count_threshold = thresholds$nc_raw_count_threshold
)

# Add flags to plate cell QC table
plate_cell_table = plate_cell_qc_flags(
  plate_cell_table = plate_cell_table,
  nc_variability_threshold = thresholds$nc_variability_threshold,
  error_rate_threshold = thresholds$error_rate_threshold,
  pc_viability_threshold = thresholds$pc_viability_threshold,
  nc_raw_count_threshold = thresholds$nc_raw_count_threshold,
  contains_poscon = contains_poscon
)

# Filter for failing cell lines
plate_cell_qc_flags_table = plate_cell_table |> dplyr::filter(!is.na(qc_flag))

# Write plate_cell_qc_table for internal use
plate_cell_int_path = file.path(args$out, "qc_tables", "plate_cell_qc_table_internal.csv")
message("Writing out internal plate_cell_qc_table to ", plate_cell_int_path)
write_out_table(table = plate_cell_table, path = plate_cell_int_path)
check_file_exists(plate_cell_int_path)

# Write plate_cell_qc_table for portal use
plate_cell_ext_path = file.path(args$out, "qc_tables", "plate_cell_qc_table.csv")
message("Writing out external plate_cell_qc_table to ", plate_cell_ext_path)
if (contains_poscon) {
  columns_to_write = c("cell_set", "pool_id", "depmap_id", "lua", "pcr_plate",
                       "pert_plate", "project_code",
                       "error_rate", "lfc_trt_poscon", "viability_trt_poscon",
                       paste0("median_raw_", negcon),
                       paste0("mad_log_normalized_", negcon),
                       paste0("median_log_normalized_", negcon),
                       paste0("n_replicates_", negcon),
                       paste0("n_replicates_", poscon),
                       "qc_pass", "qc_pass_pert_plate", "qc_flag")
} else {
  columns_to_write = c("cell_set", "pool_id", "depmap_id", "lua", "pcr_plate",
                       "pert_plate", "project_code",
                       paste0("median_raw_", negcon),
                       paste0("mad_log_normalized_", negcon),
                       paste0("median_log_normalized_", negcon),
                       paste0("n_replicates_", negcon),
                       "qc_pass", "qc_pass_pert_plate", "qc_flag")
}

# Note any expected columns that were not found - this may flag the optional columns.
missing_cols = setdiff(columns_to_write, colnames(plate_cell_table))
if (length(missing_cols) > 0) {
  message("The following columns are missing from plate_cell_qc_table: ")
  message("  ", paste(missing_cols, sep = ", "))
  warning("plate_cell_qc_table may contain missing columns!")
}
write_out_table(table = plate_cell_table |> dplyr::select(any_of(columns_to_write)),
                path = plate_cell_ext_path)
check_file_exists(plate_cell_ext_path)

# Write plate_cell_qc_flags table
plate_cell_flags_outpath = file.path(args$out, "qc_tables", "plate_cell_qc_flags.csv")
message("Writing out plate_cell_qc_flags to ", plate_cell_flags_outpath)
write_out_table(table = plate_cell_qc_flags_table, path = plate_cell_flags_outpath)
check_file_exists(plate_cell_flags_outpath)

# Create and write out cell line pert plate pass rate for portal ----
message("Creating pert plate pass rate file")
pertplate_cell_pass_rate = plate_cell_table |>
  dplyr::distinct(across(all_of(c(cell_line_cols, pert_plate_col, "qc_pass_pert_plate"))),
                  across(any_of(c("project_code")))) |>
  dplyr::group_by(across(all_of(pert_plate_col)), across(any_of(c("project_code")))) |>
  dplyr::summarise(num_cl_pass = sum(qc_pass_pert_plate),
                   num_cl_failed = sum(!qc_pass_pert_plate), .groups = "drop")
pass_rate_outpath = file.path(args$out, "qc_tables", "pertplate_cell_pass_rate.csv")
message("Writing out pertplate_cell_pass_rate to ", pass_rate_outpath)
write_out_table(pertplate_cell_pass_rate, pass_rate_outpath)
check_file_exists(pass_rate_outpath)

# End ----
message("QC module completed.")