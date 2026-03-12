options(cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(argparse)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(jsonlite)
  library(dplyr)
})
source("qc_tables/qc_tables_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser = ArgumentParser()
parser$add_argument("--normalized_counts", default = "normalized_counts.csv", help = "normalized counts file")
parser$add_argument("--sample_meta", default = "sample_meta.csv", help = "Sample meta file")
parser$add_argument("--qc_params", default = "qc_params.json", help = "File containing QC parameters")
parser$add_argument("--cell_line_cols", default = "pool_id,depmap_id,lua")
parser$add_argument("--sig_cols", default = "")
parser$add_argument("--pseudocount", default = 20)
parser$add_argument("-n", "--negcon_type", default = "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default = "trt_poscon")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")
args = parser$parse_args()

# Read in files and set up parameters ----
message("Reading in normalized_counts from", args$normalized_counts, "...")
normalized_counts = read_data_table(args$normalized_counts)
message("Reading in sample_meta from ", args$sample_meta, "...")
sample_meta = read_data_table(args$sample_meta)
message("Reading in qc_thresholds from ", args$qc_params, "...")
thresholds = load_thresholds_from_json(args$qc_params)


# Set up parameters
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
sig_cols = unlist(strsplit(args$sig_cols, ","))
pseudocount = as.numeric(args$pseudocount)
negcon = args$negcon_type
poscon = args$poscon_type
contains_poscon = any(sample_meta$pert_type == args$poscon_type)

# Create a qc output directory if one does not exist
if (!dir.exists(file.path(args$out, "qc_tables"))) {
  dir.create(file.path(args$out, "qc_tables"))
}

# Calculate plate cell QCs ----
normalized_counts_rm_cbc = filter_control_barcodes(normalized_counts)

plate_cell_table = generate_plate_cell_table(
  normalized_counts = normalized_counts_rm_cbc,
  sample_meta = sample_meta,
  cell_line_cols = cell_line_cols,
  sig_cols = sig_cols,
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
                       "error_rate", "lfc_trt_poscon",
                       "median_raw_ctl_vehicle", "mad_log_normalized_ctl_vehicle",
                       "median_log_normalized_ctl_vehicle",
                       "n_replicates_ctl_vehicle", "n_replicates_trt_poscon",
                       "viability_trt_poscon", "qc_pass", "qc_pass_pert_plate", "qc_flag")
} else {
  columns_to_write = c("cell_set", "pool_id", "depmap_id", "lua", "pcr_plate",
                       "pert_plate", "project_code",
                       "median_raw_ctl_vehicle", "mad_log_normalized_ctl_vehicle",
                       "median_log_normalized_ctl_vehicle", "n_replicates_ctl_vehicle",
                       "qc_pass", "qc_pass_pert_plate", "qc_flag")
}
# Foresee errors here with different/custom runs - wrap in try?
write_out_table(table = plate_cell_table |> dplyr::select(all_of(columns_to_write)),
                path = plate_cell_ext_path)
check_file_exists(plate_cell_ext_path)

# Write plate_cell_qc_flags table
plate_cell_flags_outpath = file.path(args$out, "qc_tables", "plate_cell_qc_flags.csv")
message("Writing out plate_cell_qc_flags to ", plate_cell_flags_outpath)
write_out_table(table = plate_cell_qc_flags_table, path = plate_cell_flags_outpath)
check_file_exists(plate_cell_flags_outpath)

# Create and write out cell line pert plate pass rate for portal ----
pass_rate_outpath = file.path(args$out, "qc_tables", "pertplate_cell_pass_rate.csv")
pertplate_cell_pass_rate = plate_cell_table |>
  dplyr::distinct(pert_plate, project_code, qc_pass_pert_plate, depmap_id, lua, cell_set) |>
  dplyr::group_by(pert_plate, project_code) |>
  dplyr::summarise(num_cl_pass = sum(qc_pass_pert_plate),
                   num_cl_failed = sum(!qc_pass_pert_plate), .groups = "drop")
pass_rate_outpath = file.path(args$out, "qc_tables", "pertplate_cell_pass_rate.csv")
message("Writing out pertplate_cell_pass_rate to ", pass_rate_outpath)
write_out_table(pertplate_cell_pass_rate, pass_rate_outpath)
check_file_exists(pass_rate_outpath)

# Calculate PCR plate QC flags ----
pcr_plate_qc_flags_table = generate_pcr_plate_qc_flags_table(
  plate_cell_table = plate_cell_table,
  fraction_expected_controls = thresholds$fraction_expected_controls,
  contains_poscon = contains_poscon
)

# Write pcr_plate_qc_flags table
pcr_plate_qc_flags_outpath = file.path(args$out, "qc_tables", "pcr_plate_qc_flags.csv")
message("Writing out pcr_plate_qc_flags to ", pcr_plate_qc_flags_outpath)
write_out_table(table = pcr_plate_qc_flags_table, path = pcr_plate_qc_flags_outpath)
check_file_exists(pcr_plate_qc_flags_outpath)

# Future filtering ----
# Filter plate_cell_qc flags from normalized counts
# plate_cell_filt_norm_counts =
#  dplyr::anti_join(normalized_counts, plate_cell_qc_flags_table,
#                   by = c("pcr_plate", "depmap_id", "lua", "pool_id"))

# Filter pcr_plate qcs from normalized counts
# final = dplyr::anti_join(plate_cell_filt_norm_counts, pcr_plate_qc_flags_table,
#                          by = c("pcr_plate", "pert_plate"))

# End ----
message("QC module completed.")