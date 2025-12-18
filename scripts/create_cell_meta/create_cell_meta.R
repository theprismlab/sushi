#!/usr/bin/env Rscript

#' Cell Metadata Creation Script
#' 
#' This script creates cell metadata files by retrieving information from CellDB
#' and processing it for use in PRISM analysis pipelines.
#'
#' @author PRISM
#' @version 1.0

# Set options and load libraries ----
options(cli.unicode = FALSE)

# Suppress package startup messages
suppressPackageStartupMessages({
  library(argparse)
  library(magrittr)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(httr)
  library(jsonlite)
})

# Load utility functions
source("utils/kitchen_utensils.R")
source("create_cell_meta/create_cell_meta_functions.R")

# Parse command line arguments ----
parser <- ArgumentParser(description = "Create cell metadata files from CellDB")

# Add arguments
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                   help = "Print detailed output [default]")
parser$add_argument("-q", "--quietly", action = "store_false",
                   dest = "verbose", help = "Print minimal output")
parser$add_argument("--wkdir", default = getwd(), 
                   help = "Working directory")
parser$add_argument("-o", "--out", default = "", 
                   help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default = "sample_meta.csv", 
                   help = "Sample metadata file path")
parser$add_argument("--cb_ladder", default = "", 
                   help = "Control barcode ladder name or file path")
parser$add_argument("--api_url", default = "https://api.clue.io/api/", 
                   help = "Base API URL for CellDB")
parser$add_argument("--api_key", default = "", 
                   help = "Clue API key")
parser$add_argument("--cell_line_cols", default = "depmap_id,lua,pool_id,depmap_id,lua,cell_set,growth_condition",
                   help = "Comma-separated list of cell line metadata columns to include")

# Parse arguments
args <- parser$parse_args()

# Set cell_line_cols
cell_line_cols <- unlist(strsplit(args$cell_line_cols, split = ","))

# Setup output directory ----
if (args$out == "") {
  args$out <- args$wkdir
}

if (args$verbose) {
  message("=== Cell Metadata Creation ===")
  message(sprintf("Working directory: %s", args$wkdir))
  message(sprintf("Output directory: %s", args$out))
}

# Validate control barcode ladder ----
if (args$cb_ladder == "") {
  stop("ERROR: Control barcode ladder not specified in arguments.")
} else {
  cb_ladder <- args$cb_ladder
  if (args$verbose) {
    message(sprintf("Using control barcode ladder: %s", cb_ladder))
  }
}

# Load sample metadata ----
if (!file.exists(args$sample_meta)) {
  stop(sprintf("ERROR: Sample metadata file not found: %s", args$sample_meta))
}

sample_meta <- read.csv(args$sample_meta)

if (args$verbose) {
  message(sprintf("Loaded sample metadata with %d rows", nrow(sample_meta)))
}

# Get API key ----
if (args$api_key != "") {
  api_key <- args$api_key
} else if (Sys.getenv("API_KEY") != "") {
  api_key <- Sys.getenv("API_KEY")
  if (args$verbose) {
    message("Using API key from environment variable")
  }
} else {
  stop("ERROR: No API key provided via argument or environment variable API_KEY.")
}

# Fetch data from CellDB ----
if (args$verbose) {
  message("Fetching cell information from CellDB...")
  message(sprintf("API URL: %s", args$api_url))
}

# Construct API URLs
cell_sets_url <- paste(args$api_url, "cell-db", "cell-sets", sep = "/")
cell_pools_url <- paste(args$api_url, "cell-db", "assay-pools", sep = "/")
cell_lines_url <- paste(args$api_url, "cell-db", "cell-lines", sep = "/")
assay_pools_url <- paste(args$api_url, "cell-db", "cell-sets", sep = "/")
control_barcodes_url <- paste(args$api_url, "cell-db", "bq-control-barcodes", sep = "/")

# Fetch data
cell_sets_df <- get_cell_api_info(cell_sets_url, api_key)
cell_pools_df <- get_cell_api_info(cell_pools_url, api_key)
cell_lines_df <- get_cell_api_info(cell_lines_url, api_key)
assay_pools_df <- get_cell_api_info(assay_pools_url, api_key)
control_barcodes_df <- get_cell_api_info(control_barcodes_url, api_key)

if (is.null(cell_sets_df) || is.null(cell_pools_df) || 
    is.null(cell_lines_df) || is.null(assay_pools_df)) {
  stop("ERROR: Failed to retrieve data from CellDB API")
}

if (args$verbose) {
  message(sprintf("Retrieved %d rows from cell-sets", nrow(cell_sets_df)))
  message(sprintf("Retrieved %d rows from cell-pools", nrow(cell_pools_df)))
  message(sprintf("Retrieved %d row from cell-lines", nrow(cell_lines_df)))
}

# Process assay pools metadata
assay_pools_meta <- dplyr::select(assay_pools_df, -cell_set_desc)

# Get control barcode data ----
if (args$verbose) {
  message("Retrieving control barcode data...")
}

if (!str_detect(cb_ladder, ".csv")) {
  # Fetch control barcode data from API
  if (args$verbose) {
    message(sprintf("Fetching control barcode data for ladder: %s", cb_ladder))
  }
  
  control_bc_url <- paste(args$api_url, "v_control_barcodes", sep = "/")
  control_bc_df <- get_cell_api_info(
    control_bc_url,
    api_key,
    filter = list(
      where = list(set = cb_ladder)
    )
  )
  
  if (is.null(control_bc_df) || nrow(control_bc_df) == 0) {
    stop(sprintf("ERROR: No control barcode data found for ladder: %s", cb_ladder))
  }
  
  # Rename and select columns
  control_bc_df <- control_bc_df %>%
    dplyr::rename(
      cb_name = name,
      forward_read_barcode = sequence,
      cb_log10_dose = log_dose,
      cb_ladder = set
    ) %>%
    dplyr::mutate(cb_ladder = tolower(cb_ladder)) %>%
    dplyr::select(forward_read_barcode, cb_name, cb_log10_dose, cb_ladder)
  
} else {
  # Load control barcode data from file
  file_path <- file.path(args$out, cb_ladder)
  
  if (!file.exists(file_path)) {
    stop(sprintf("ERROR: Control barcode file not found: %s", file_path))
  }
  
  if (args$verbose) {
    message(sprintf("Loading control barcode data from file: %s", file_path))
  }
  
  control_bc_df <- read.csv(file_path)
}

if (args$verbose) {
  message(sprintf("Retrieved %d control barcodes", nrow(control_bc_df)))
}

# Prepare cell line metadata ----
if (args$verbose) {
  message("Preparing cell line metadata...")
}

cell_line_meta_cols <- c('depmap_id', 'forward_read_barcode', 'lua')
cell_line_meta <- cell_lines_df %>%
  dplyr::rename("forward_read_barcode" = "dna_sequence") %>% 
  dplyr::select(dplyr::any_of(c(cell_line_meta_cols))) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(depmap_id)) %>%
    dplyr::filter(!is.na(forward_read_barcode))

# Set LUA to null for control barcodes
if (nrow(control_barcodes_df) > 0) {
    cell_line_meta <- cell_line_meta %>%
    dplyr::mutate(depmap_id = ifelse(forward_read_barcode %in% control_barcodes_df$dna_sequence, NA, depmap_id))
}
if (args$verbose) {
  message(sprintf("Found metadata for %d cell lines", nrow(cell_line_meta)))
}

# Create cell set metadata ----
if (args$verbose) {
  message("Creating cell set metadata...")
}

cell_sets_result <- create_cell_set_meta(
  sample_meta, 
  cell_sets_df, 
  cell_pools_df, 
  cell_line_meta,
  api_key,
  args$api_url,
  args$verbose
)

cell_set_meta <- cell_sets_result[[1]]
failed_cell_sets <- cell_sets_result[[2]]

# Check for failed cell sets
if (length(failed_cell_sets) > 0) {
  message("ERROR: The following cell sets are not registered in CellDB:")
  for (failed_set in failed_cell_sets) {
    message(sprintf("  - %s", failed_set))
  }
  stop("Cannot proceed with unregistered cell sets")
}

# Convert cell set metadata to long form ----
if (args$verbose) {
  message("Converting cell set metadata to long form...")
}

cell_set_meta_long <- cell_set_meta %>%
  tidyr::separate_rows(members, sep = ";")

if (args$verbose) {
  message(sprintf("Created long-form cell set metadata with %d total cell lines", nrow(cell_set_meta_long)))
}

# Combine with assay pool metadata ----
if (args$verbose) {
  message("Combining cell set metadata with assay pool metadata...")
}

# Check if all cell sets exist in assay pools
all_sets_exist <- all(cell_set_meta_long$cell_set %in% assay_pools_meta$davepool_id)

if (all_sets_exist) {
  if (args$verbose) {
    message("All cell sets found in assay pool metadata")
  }
  
  cell_set_assay_pool_meta <- cell_set_meta_long %>%
    dplyr::inner_join(assay_pools_meta, by = c("cell_set" = "davepool_id", "members" = "depmap_id"),
                      relationship = "many-to-many") %>%
    dplyr::select(cell_set, pool_id, barcode_id, growth_condition, depmap_id = members) %>%
    dplyr::rename("lua" = "barcode_id") %>%
    dplyr::distinct()

  # Make sure all rows contain a value for growth_condition. If not, stop with error.
  if (any(is.na(cell_set_assay_pool_meta$growth_condition)) && "growth_condition" %in% cell_line_cols) {
    stop(paste(
      "ERROR: One or more cell lines are missing growth_condition values. Missing lines are:",
      paste(cell_set_assay_pool_meta %>% filter(is.na(growth_condition)) %>% pull(depmap_id), collapse = ", "),
      "Please ensure all cell lines have growth_condition values assigned in cellDB, or remove growth_condition from cell_line_cols, disable the bias correction module and re-run.",
      sep = " "
    ))
  }
  
} else {
  message("WARNING: One or more cell sets not found in assay pool metadata")
  message("Unable to include pool_id in cell set metadata")
  
  cell_set_assay_pool_meta <- cell_set_meta_long %>%
    dplyr::select(cell_set, depmap_id = members) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(cell_line_meta %>% select(depmap_id, lua), by = "depmap_id")
}

# Prepare control barcode metadata ----
if (args$verbose) {
  message("Preparing control barcode metadata...")
}

cb_meta <- control_bc_df

# Add depmap_id and lua if missing
if (!all(c("depmap_id", "lua") %in% colnames(cb_meta))) {
  if (args$verbose) {
    message("Adding lua column to control barcode metadata")
  }

  cb_meta <- cb_meta %>%
    dplyr::left_join(cell_line_meta %>% select(c("forward_read_barcode","lua")), by = "forward_read_barcode")
}

# Add null depmap_ids
cb_meta <- cb_meta %>%
  dplyr::mutate(depmap_id = NA)

# Write output files ----
if (args$verbose) {
  message("Writing output files...")
}

# Write cell line metadata
cell_line_out_file <- file.path(args$out, 'cell_line_meta.csv')
if (args$verbose) {
  message(sprintf("Writing cell line metadata to: %s", cell_line_out_file))
}
write_out_table(cell_line_meta, cell_line_out_file)


# Write cell set and pool metadata
cell_set_assay_pool_out_file <- file.path(args$out, 'cell_set_and_pool_meta.csv')
if (args$verbose) {
  message(sprintf("Writing cell set and pool metadata to: %s", cell_set_assay_pool_out_file))
}
write_out_table(cell_set_assay_pool_meta, cell_set_assay_pool_out_file)

# Write control barcode metadata if available
if (nrow(control_bc_df) > 0) {
  control_barcode_out_file <- file.path(args$out, 'CB_meta.csv')
  if (args$verbose) {
    message(sprintf("Writing control barcode metadata to: %s", control_barcode_out_file))
  }
  write_out_table(cb_meta, control_barcode_out_file)
}

# Verify output files ----
check_file_exists(cell_line_out_file)
check_file_exists(cell_set_assay_pool_out_file)
# check_file_exists(portal_cell_line_out_file)

if (args$verbose) {
  message("=== Cell Metadata Creation Complete ===")
}
