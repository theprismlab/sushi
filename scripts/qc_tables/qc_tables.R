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
parser <- ArgumentParser()
parser$add_argument(
  "--cell_set_and_pool_meta",
  default = "cell_set_and_pool_meta.csv",
  help = "Cell line metadata"
)
parser$add_argument(
  "--normalized_counts",
  default = "normalized_counts.csv", help = "normalized counts file"
)
parser$add_argument(
  "--annotated_counts",
  default = "annotated_counts.csv", help = "annotated counts file"
)
parser$add_argument(
  "--filtered_counts",
  default = "filtered_counts.csv", help = "filtered counts file"
)
parser$add_argument(
  "-o", "--out",
  default = getwd(), help = "Output path. Default is working directory"
)
parser$add_argument(
  "--control_barcode_meta",
  default = "CB_meta.csv", help = "Control barcode metadata"
)
parser$add_argument(
  "--unknown_barcode_counts",
  default = "unknown_barcode_counts.csv",
  help = "Unknown barcode counts file"
)
parser$add_argument("-n", "--negcon_type", default = "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default = "trt_poscon")
parser$add_argument("--cell_line_cols", default = "depmap_id,pool_id,lua")
parser$add_argument("--id_cols", default = "pcr_plate,pcr_well")
parser$add_argument("--sig_cols", default = "")
parser$add_argument("--count_threshold", default = 40)
parser$add_argument("--pseudocount", default = 20)
parser$add_argument("--filter_qc_flags",
  default = "true",
  help = "Filter out wells with QC flags. Default is TRUE"
)
parser$add_argument(
  "--qc_params",
  default = "qc_params.json",
  help = "File containing QC parameters"
)
parser$add_argument(
  "--sample_meta",
  default = "sample_meta.csv",
  help = "Sample metadata file"
)
parser$add_argument(
  "--prism_barcode_counts",
  default = "prism_barcode_counts.csv",
  help = "Prism barcode counts file"
)
parser$add_argument(
  "--cell_line_meta",
    default = "cell_line_meta.csv",
    help = "Cell line metadata file containing all PRISM cell lines"
)

args <- parser$parse_args()

# Read in files as data.table objects ----
print(paste0("Reading in ", args$cell_set_and_pool_meta, "....."))
cell_set_meta <- read_data_table(args$cell_set_and_pool_meta)
print(paste0("Reading in ", args$annotated_counts, "....."))
annotated_counts <- read_data_table(args$annotated_counts)
print(paste0("Reading in ", args$filtered_counts, "....."))
filtered_counts <- read_data_table(args$filtered_counts)
print(paste0("Reading in ", args$control_barcode_meta, "....."))
cb_meta <- read_data_table(args$control_barcode_meta)
print(paste0("Reading in ", args$unknown_barcode_counts, "....."))
unknown_counts <- read_data_table(args$unknown_barcode_counts)
print(paste0("Reading in ", args$sample_meta, "....."))
sample_meta <- read_data_table(args$sample_meta)
# If normzlied_counts_original.csv exists, use that, otherwise use args$normalized_counts
normalized_counts_original_path <- paste0(args$out, "/normalized_counts_original.csv")
if (file.exists(normalized_counts_original_path)) {
  print("Original normalized counts found, will use this file for QC flags and tables.")
  print("Reading in normalized_counts_original.csv.....")
  normalized_counts <- read_data_table(normalized_counts_original_path)
} else {
print(paste0("Reading in ", args$normalized_counts, "....."))
normalized_counts <- read_data_table(args$normalized_counts)
}
print(paste0("Reading in ", args$prism_barcode_counts, "....."))
prism_barcode_counts <- read_data_table(args$prism_barcode_counts)
print(paste0("Reading in ", args$cell_line_meta, "....."))
cell_line_meta <- read_data_table(args$cell_line_meta)

# Check if the output directory exists, if not create it
if (!dir.exists(paste0(args$out, "/qc_tables"))) {
  dir.create(paste0(args$out, "/qc_tables"))
}

# Check if any samples are positive controls, if not set poscon to FALSE, if so TRUE
contains_poscon <- any(sample_meta$pert_type == args$poscon_type)

# DEFINE COLUMNS
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
id_cols <- args$id_cols
id_cols_list <- strsplit(id_cols, ",")[[1]]
sig_cols = unlist(strsplit(args$sig_cols, ","))
count_threshold <- as.numeric(args$count_threshold)
pseudocount <- as.numeric(args$pseudocount)
filter_qc_flags <- as.logical(toupper(args$filter_qc_flags))

# Set control types
poscon <- args$poscon_type
negcon <- args$negcon_type

# LOAD QC PARAMETERS
thresholds <- load_thresholds_from_json(args$qc_params)


# Calculate the number of expected poscons and negcons
n_expected_controls <- sample_meta %>%
  dplyr::filter(pert_type %in% c(negcon, poscon)) %>%
  group_by(pert_plate, pcr_plate, pert_type) %>%
  summarize(unique_bio_rep = n_distinct(bio_rep), .groups = "drop") %>%
  pivot_wider(
    names_from = pert_type,
    values_from = unique_bio_rep,
    names_prefix = "n_expected_"
  )

# ID COLS
id_cols_table <- generate_id_cols_table(
  normalized_counts = normalized_counts, annotated_counts = annotated_counts, unknown_counts = unknown_counts,
  cell_set_meta = cell_set_meta, id_cols_list = id_cols_list, cell_line_cols = cell_line_cols,
  count_threshold = count_threshold, cb_meta = cb_meta, pseudocount = pseudocount
)

id_cols_qc_flags_table <- id_cols_qc_flags(id_cols_table = id_cols_table,
                                           contamination_threshold = thresholds$contamination_threshold,
                                           cb_mae_threshold = thresholds$cb_mae_threshold,
                                           cb_spearman_threshold = thresholds$cb_spearman_threshold,
                                           cb_cl_ratio_low_negcon = thresholds$cb_cl_ratio_low_negcon,
                                           cb_cl_ratio_high_poscon = thresholds$cb_cl_ratio_high_poscon,
                                           cb_cl_ratio_low_poscon = thresholds$cb_cl_ratio_low_poscon,
                                           cb_cl_ratio_high_negcon = thresholds$cb_cl_ratio_high_negcon,
                                           well_reads_threshold = thresholds$well_reads_threshold)

id_cols_filtered_normalized_counts <- dplyr::anti_join(normalized_counts, id_cols_qc_flags_table, by = c("pcr_plate", "pcr_well"))


# POOL WELL
pool_well_table <- generate_pool_well_qc_table(
  normalized_counts = id_cols_filtered_normalized_counts,
  pool_well_delta_threshold = thresholds$pool_well_delta_threshold,
  pool_well_fraction_threshold = thresholds$pool_well_fraction_threshold
)

pool_well_qc_flags_table <- pool_well_table %>%
  dplyr::filter(!is.na(qc_flag)) %>%
  unique()

pool_well_filtered_normalized_counts <- dplyr::anti_join(id_cols_filtered_normalized_counts, pool_well_qc_flags_table, by = c("pcr_plate", "pcr_well", "pool_id"))


# PLATE CELL

filtered_normalized_counts_rm_ctl <- pool_well_filtered_normalized_counts %>%
  filter_control_barcodes()
filtered_counts_rm_ctl <- filtered_counts %>%
  filter_control_barcodes()

plate_cell_table <- generate_cell_plate_table(
  normalized_counts = filtered_normalized_counts_rm_ctl, filtered_counts = filtered_counts_rm_ctl,
  cell_line_cols = cell_line_cols, sig_cols = sig_cols, pseudocount = pseudocount, contains_poscon = contains_poscon,
  poscon = poscon, negcon = negcon,
  nc_variability_threshold = thresholds$nc_variability_threshold,
  error_rate_threshold = thresholds$error_rate_threshold,
  pc_viability_threshold = thresholds$pc_viability_threshold,
  nc_raw_count_threshold = thresholds$nc_raw_count_threshold
)

plate_cell_qc_flags_table <- plate_cell_qc_flags(
  plate_cell_table = plate_cell_table,
  nc_variability_threshold = thresholds$nc_variability_threshold,
  error_rate_threshold = thresholds$error_rate_threshold,
  pc_viability_threshold = thresholds$pc_viability_threshold,
  nc_raw_count_threshold = thresholds$nc_raw_count_threshold,
  contains_poscon = contains_poscon) %>%
  dplyr::filter(!is.na(qc_flag))


plate_cell_filtered_normalized_counts <-
  dplyr::anti_join(
  pool_well_filtered_normalized_counts, plate_cell_qc_flags_table,
  by = c("pcr_plate", "depmap_id", "lua", "pool_id")
)


# PCR PLATE

pcr_plate_qc_flags_table <- generate_pcr_plate_qc_flags_table(plate_cell_table = plate_cell_table,
                                                              fraction_expected_controls = thresholds$fraction_expected_controls,
                                                              contains_poscon = contains_poscon)

final_filtered_normalized_counts <- pool_well_filtered_normalized_counts
#  dplyr::anti_join(
#    plate_cell_filtered_normalized_counts, pcr_plate_qc_flags_table,
#    by = c("pcr_plate", "pert_plate"))

# WRITE OUT RESULTS --------

# Write pcr_plate_qc_flags table
pcr_plate_qc_flags_outpath <- paste0(args$out, "/qc_tables/pcr_plate_qc_flags.csv")
print(paste0("Writing out pcr_plate_qc_flags to ", pcr_plate_qc_flags_outpath))
write_out_table(table = pcr_plate_qc_flags_table, path = pcr_plate_qc_flags_outpath)
check_file_exists(pcr_plate_qc_flags_outpath)

# Write plate_cell_qc_table for internal use
plate_cell_qc_table_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table_internal.csv")
print(paste0("Writing out internal plate_cell_qc_table to ", plate_cell_qc_table_outpath))
write_out_table(table = plate_cell_table, path = plate_cell_qc_table_outpath)
check_file_exists(plate_cell_qc_table_outpath)

# Write plate_cell_qc_table for portal use
plate_cell_qc_table_outpath_external <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out external plate_cell_qc_table to ", plate_cell_qc_table_outpath_external))
if (contains_poscon) {
  columns_to_write <- c(
    "cell_set","pool_id", "depmap_id", "lua", "pcr_plate",
    "pert_plate", "project_code",
    "error_rate", "lfc_trt_poscon",
    "median_raw_ctl_vehicle", "mad_log_normalized_ctl_vehicle",
    "median_log_normalized_ctl_vehicle",
    "n_replicates_ctl_vehicle", "n_replicates_trt_poscon",
    "viability_trt_poscon", "qc_pass", "qc_pass_pert_plate"
  )
  } else {
    columns_to_write <- c(
    "cell_set","pool_id", "depmap_id", "lua", "pcr_plate",
    "pert_plate", "project_code",
    "median_raw_ctl_vehicle", "mad_log_normalized_ctl_vehicle",
    "median_log_normalized_ctl_vehicle",
    "n_replicates_ctl_vehicle",
    "qc_pass", "qc_pass_pert_plate"
  )
}
write_out_table(table = plate_cell_table %>% dplyr::select(all_of(columns_to_write)),
                path = plate_cell_qc_table_outpath_external)
check_file_exists(plate_cell_qc_table_outpath_external)

# Write cl_pertplate_pass_rate_qc_table for portal
pertplate_cell_pass_rate_outpath <- paste0(args$out, "/qc_tables/pertplate_cell_pass_rate.csv")
pertplate_cell_pass_rate <- plate_cell_table %>% 
  dplyr::distinct(pert_plate, project_code, qc_pass_pert_plate, depmap_id, lua, cell_set) %>%
  dplyr::group_by(pert_plate, project_code) %>%
  dplyr::summarise(num_cl_pass=sum(qc_pass_pert_plate),
                   num_cl_failed=sum(!qc_pass_pert_plate)) %>%
  dplyr::ungroup()
write_out_table(pertplate_cell_pass_rate, pertplate_cell_pass_rate_outpath)
check_file_exists(pertplate_cell_pass_rate_outpath)

# Write plate_cell_qc_flags table
plate_cell_qc_flags_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_flags.csv")
print(paste0("Writing out plate_cell_qc_flags to ", plate_cell_qc_flags_outpath))
write_out_table(table = plate_cell_qc_flags_table, path = plate_cell_qc_flags_outpath)
check_file_exists(plate_cell_qc_flags_outpath)

# Write pool_well_qc_table ----------
pool_well_qc_table_outpath <- paste0(args$out, "/qc_tables/pool_well_qc_table.csv")
print(paste0("Writing out pool_well_qc_table to ", pool_well_qc_table_outpath))
write_out_table(table = pool_well_table, path = pool_well_qc_table_outpath)
check_file_exists(pool_well_qc_table_outpath)

# Write pool_well_qc_flags table ----------
pool_well_qc_flags_outpath <- paste0(args$out, "/qc_tables/pool_well_qc_flags.csv")
print(paste0("Writing out pool_well_qc_flags to ", pool_well_qc_flags_outpath))
write_out_table(table = pool_well_qc_flags_table, path = pool_well_qc_flags_outpath)
check_file_exists(pool_well_qc_flags_outpath)

# Write id_cols table ----------
id_cols_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_table.csv")
print(paste0("Writing out id_cols_qc_table to ", id_cols_outpath))
write_out_table(table = id_cols_table, path = id_cols_outpath)
check_file_exists(id_cols_outpath)

# Write id_cols_qc_flags table ----------
id_cols_qc_flags_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_flags.csv")
print(paste0("Writing out id_cols_qc_flags to ", id_cols_qc_flags_outpath))
write_out_table(table = id_cols_qc_flags_table, path = id_cols_qc_flags_outpath)
check_file_exists(id_cols_qc_flags_outpath)

if (args$filter_qc_flags) {
  # Filter out wells with QC flags
  print("Filtering out wells with QC flags")
  # Write original normalized counts ----------
  normalized_counts_original_outpath <- paste0(args$out, "/normalized_counts_original.csv")
  print(paste0("Writing unfiltered normalized_counts to ", normalized_counts_original_outpath))
  write_out_table(table = normalized_counts, path = normalized_counts_original_outpath)
  check_file_exists(normalized_counts_original_outpath)

  # Write filtered normalized counts ----------
  filtered_normalized_counts_outpath <- paste0(args$out, "/normalized_counts.csv")
  print(paste0("Writing filtered normalized_counts to ", filtered_normalized_counts_outpath))
  write_out_table(table = final_filtered_normalized_counts, path = filtered_normalized_counts_outpath)
  check_file_exists(filtered_normalized_counts_outpath)
} else {
  print("Nomalized counts not filtered for qc_flags.")
}

# Compute variance decomposition table
print("Computing variance decomposition table...")
variance_decomp <- compute_variance_decomposition(
  normalized_counts = normalized_counts,
  metric = "n",
  cell_line_cols = cell_line_cols,
  id_cols = id_cols_list,
  negcon = negcon
)
variance_decomp_outpath <- paste0(args$out, "/qc_tables/variance_decomposition.csv")
print(paste0("Writing out variance_decomposition to ", variance_decomp_outpath))
write_out_table(table = variance_decomp, path = variance_decomp_outpath)
check_file_exists(variance_decomp_outpath)

# Compute contamination tables
contamination_tables <- compute_contamination_qc_tables(
  prism_barcode_counts = prism_barcode_counts,
  unknown_barcode_counts = unknown_counts,
  cell_set_and_pool_meta = cell_set_meta,
  cell_line_meta = cell_line_meta,
  cb_meta = cb_meta,
  sample_meta = sample_meta
)

# Write out contamination tables
for (table_name in names(contamination_tables)) {
  table_outpath <- paste0(args$out, "/qc_tables/", table_name, ".csv")
  print(paste0("Writing out ", table_name, " to ", table_outpath))
  write_out_table(table = contamination_tables[[table_name]], path = table_outpath)
  check_file_exists(table_outpath)
}

paste0("QC module completed.")
