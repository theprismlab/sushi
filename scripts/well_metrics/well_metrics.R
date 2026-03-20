# Script for well level QC filters

options(cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(data.table)
  library(tidyverse)
})

source("well_metrics/well_metrics_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser = ArgumentParser()

# Upstream Sushi outputs
parser$add_argument("--unknown_barcode_counts", default = "unknown_barcode_counts.csv",
                    help = "Unknown barcode counts file")
parser$add_argument("--prism_barcode_counts", default = "prism_barcode_counts.csv",
                    help = "PRISM barcode counts file")
parser$add_argument("--annotated_counts", default = "annotated_counts.csv", help = "annotated counts file")
parser$add_argument("--normalized_counts", default = "normalized_counts.csv", help = "normalized counts file")

# Sushi metadata
parser$add_argument("--cell_set_and_pool_meta", default = "cell_set_and_pool_meta.csv", help = "Cell line metadata")
parser$add_argument("--cb_meta", default = "CB_meta.csv", help = "Control barcode metadata")
parser$add_argument("--cell_line_meta", default = "cell_line_meta.csv", help = "Cell line metadata")
parser$add_argument("--sample_meta", default = "sample_meta.csv", help = "Sample metadata")
parser$add_argument("--qc_params", default = "qc_params.json", help = "File containing QC parameters")

# Other params
parser$add_argument("--id_cols", default = "pcr_plate,pcr_well")
parser$add_argument("--cell_line_cols", default = "depmap_id,pool_id,lua")
parser$add_argument("--negcon_type", default = "ctl_vehicle")
parser$add_argument("--count_threshold", default = 40)
parser$add_argument("--pseudocount", default = 20)
parser$add_argument("--filter_qc_flags", default = "true", help = "Filter out wells with QC flags. Default is TRUE")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")

args = parser$parse_args()

# Read in files and set up parameters ----
message("Reading in unknown_barcode_counts from ", args$unknown_barcode_counts, "...")
unknown_counts = read_data_table(args$unknown_barcode_counts)
message("Reading in prism_barcode_counts from ", args$prism_barcode_counts, "...")
prism_barcode_counts = read_data_table(args$prism_barcode_counts)
message("Reading in annotated_counts from ", args$annotated_counts, "...")
annotated_counts = read_data_table(args$annotated_counts)

# Read in normalized_counts_original if possible
norm_counts_original_path = file.path(args$out, "normalized_counts_original.csv")
if (file.exists(norm_counts_original_path)) {
  message("Found normalized_counts_original.csv. This file is preferred for well level metrics and QC flags.")
  message("Reading in normalized_counts_original from ", norm_counts_original_path, "...")
  normalized_counts = read_data_table(norm_counts_original_path)
} else {
  message("Reading in normalized_counts from", args$normalized_counts, "...")
  normalized_counts = read_data_table(args$normalized_counts)
}

message("Reading in cell_set_and_pool_meta from ", args$cell_set_and_pool_meta, "...")
cell_set_and_pool_meta = read_data_table(args$cell_set_and_pool_meta)
message("Reading in cb_meta from ", args$cb_meta, "...")
cb_meta = read_data_table(args$cb_meta)
message("Reading in cell_line_meta from ", args$cell_line_meta, "...")
cell_line_meta = read_data_table(args$cell_line_meta)
message("Reading in sample_meta from ", args$sample_meta, "...")
sample_meta = read_data_table(args$sample_meta)
message("Reading in qc_thresholds from ", args$qc_params, "...")
thresholds = load_thresholds_from_json(args$qc_params)

# Set up parameters
id_cols = unlist(strsplit(args$id_cols, ","))
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
count_threshold = as.numeric(args$count_threshold)
pseudocount = as.numeric(args$pseudocount)
filter_qc_flags = as.logical(toupper(args$filter_qc_flags))
negcon = args$negcon_type

# Create a qc output directory if one does not exist
if (!dir.exists(file.path(args$out, "qc_tables"))) {
  dir.create(file.path(args$out, "qc_tables"))
}

# ID COLS QC ----
message("Calculating id_cols level QCs ...")
id_cols_table = generate_id_cols_table(
  unknown_counts = unknown_counts,
  annotated_counts = annotated_counts,
  normalized_counts = normalized_counts,
  cb_meta = cb_meta,
  id_cols = id_cols,
  cell_line_cols = cell_line_cols,
  pseudocount = pseudocount,
  count_threshold = count_threshold,
  cb_threshold = thresholds$well_reads_threshold,
  cb_spearman_threshold = thresholds$cb_spearman_threshold,
  cb_mae_threshold = thresholds$cb_mae_threshold
)

# Identify poor performing wells to drop
id_cols_qc_flags_table = id_cols_qc_flags(
  id_cols_table = id_cols_table,
  contamination_threshold = thresholds$contamination_threshold,
  cb_mae_threshold = thresholds$cb_mae_threshold,
  cb_spearman_threshold = thresholds$cb_spearman_threshold,
  cb_cl_ratio_low_negcon = thresholds$cb_cl_ratio_low_negcon,
  cb_cl_ratio_high_poscon = thresholds$cb_cl_ratio_high_poscon,
  cb_cl_ratio_low_poscon = thresholds$cb_cl_ratio_low_poscon,
  cb_cl_ratio_high_negcon = thresholds$cb_cl_ratio_high_negcon,
  well_reads_threshold = thresholds$well_reads_threshold
)

# Filter normalize counts
id_cols_filt_normalized_counts = dplyr::anti_join(normalized_counts, id_cols_qc_flags_table, by = id_cols)

# Write out id col qc and flags
id_cols_outpath = file.path(args$out, "qc_tables", "id_cols_qc_table.csv")
message("Writing out id_cols_qc_table to ", id_cols_outpath)
write_out_table(table = id_cols_table, path = id_cols_outpath)
check_file_exists(id_cols_outpath)

id_cols_qc_flags_outpath = file.path(args$out, "qc_tables", "id_cols_qc_flags.csv")
message("Writing out id_cols_qc_flags to ", id_cols_qc_flags_outpath)
write_out_table(table = id_cols_qc_flags_table, path = id_cols_qc_flags_outpath)
check_file_exists(id_cols_qc_flags_outpath)

# POOL WELL QCs ----
message("Calculating pool_well level QCs ...")
pool_well_table = generate_pool_well_qc_table(
  normalized_counts = id_cols_filt_normalized_counts,
  pool_well_delta_threshold = thresholds$pool_well_delta_threshold,
  pool_well_fraction_threshold = thresholds$pool_well_fraction_threshold
)

# Filter out poor performing pools in wells
pool_well_qc_flags_table = pool_well_table |>
  dplyr::filter(!is.na(qc_flag)) |>
  unique()

pool_well_filt_norm_counts = dplyr::anti_join(id_cols_filt_normalized_counts, pool_well_qc_flags_table,
                                              by = c(id_cols, "pool_id"))

# Write out pool_well qcs and flags
pool_well_qc_table_outpath = file.path(args$out, "qc_tables", "pool_well_qc_table.csv")
message("Writing out pool_well_qc_table to ", pool_well_qc_table_outpath)
write_out_table(table = pool_well_table, path = pool_well_qc_table_outpath)
check_file_exists(pool_well_qc_table_outpath)

pool_well_qc_flags_outpath = file.path(args$out, "qc_tables", "pool_well_qc_flags.csv")
message("Writing out pool_well_qc_flags to ", pool_well_qc_flags_outpath)
write_out_table(table = pool_well_qc_flags_table, path = pool_well_qc_flags_outpath)
check_file_exists(pool_well_qc_flags_outpath)

# Write out filtered normalized counts ----
if (filter_qc_flags == TRUE) {
  # Filter out wells with QC flags
  message("Filtering out wells with QC flags from normalized counts ...")
  norm_counts_original_outpath = file.path(args$out, "normalized_counts_original.csv")
  message("Writing unfiltered normalized_counts to ", norm_counts_original_outpath)
  write_out_table(table = normalized_counts, path = norm_counts_original_outpath)
  check_file_exists(norm_counts_original_outpath)

  filt_norm_counts_outpath = file.path(args$out, "normalized_counts.csv")
  message("Writing filtered normalized_counts to ", filt_norm_counts_outpath)
  write_out_table(table = pool_well_filt_norm_counts, path = filt_norm_counts_outpath)
  check_file_exists(filt_norm_counts_outpath)
} else {
  message("Normalized counts NOT filtered for qc_flags.")
}

# VARIANCE DECOMPOSITION ----
print("Computing variance decomposition table...")
variance_decomp <- compute_variance_decomposition(
  normalized_counts = normalized_counts,
  id_cols = id_cols,
  cell_line_cols = cell_line_cols,
  metric = "n",
  negcon =  negcon
)
variance_decomp_outpath <- file.path(args$out, "qc_tables", "variance_decomposition.csv")
message("Writing out variance_decomposition to ", variance_decomp_outpath)
write_out_table(table = variance_decomp, path = variance_decomp_outpath)
check_file_exists(variance_decomp_outpath)

# CONTAMINATION TABLES ----
contamination_tables <- compute_contamination_qc_tables(
  prism_barcode_counts = prism_barcode_counts,
  unknown_barcode_counts = unknown_counts,
  cell_set_and_pool_meta = cell_set_and_pool_meta,
  cell_line_meta = cell_line_meta,
  cb_meta = cb_meta,
  sample_meta = sample_meta,
  id_cols = id_cols
)

# Write out contamination tables
for (table_name in names(contamination_tables)) {
  table_outpath <- file.path(args$out, "qc_tables", paste0(table_name, ".csv"))
  message("Writing out ", table_name, " to ", table_outpath)
  write_out_table(table = contamination_tables[[table_name]], path = table_outpath)
  check_file_exists(table_outpath)
}

message("finished well metric calculations")