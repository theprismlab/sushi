options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("biomarker/biomarker_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# Data paths
parser$add_argument("--biomarker_file", default="/data/biomarker/current/prism_biomarker_public_0625.h5",
                    help="File containing depmap data")
parser$add_argument("--collapsed_lfc", default = "collapsed_l2fc.csv",
                    help = "file containing replicate collapsed log fold change values")
parser$add_argument("--drc_file", default="DRC_TABLE.csv", help="File containing auc values from dose response")
parser$add_argument("--synergy_path", default = "synergy_scores.csv",
                    help = "File containing auc values from dose response")

# Columns and other params
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id,lua",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--collapsed_l2fc_column", default="median_l2fc",
                    help = "column containing replicate collapsed log fold change values")
parser$add_argument("--dr_column", default="log2_auc",
                    help = "column containing auc values from dose response")
parser$add_argument("--synergy_col", default = "synergy",
                    help = "Column in the synergy table containing the synergy scores.")
parser$add_argument("--univariate_biomarker", default="true", help="Whether to calculate univariate biomarkers")
parser$add_argument("--multivariate_biomarker", default="true", help="Whether to calculate multivariate biomarkers")
parser$add_argument("--lfc_biomarker", default="true", help="Whether to calculate lfc biomarkers")
parser$add_argument("--auc_biomarker", default="true", help="Whether to calculate auc biomarkers")
parser$add_argument("--synergy_biomarker", default = "false",
                    help = "Toggle to calculate biomarkers using synergy scores.")

# Paths
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--out_path", default= "", help = "Path to the output directory (per comppound)")

#TODO: Update column headers

# Get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Basic parameters ----
sig_cols = unlist(strsplit(args$sig_cols, ","))
build_dir = args$build_dir
bio_file = args$biomarker_file

# Paths and columns to use
lfc_path = args$collapsed_lfc
lfc_column = args$collapsed_l2fc_column
dr_column = args$dr_column
drc_file = args$drc_file
synergy_col = args$synergy_col
synergy_path = args$synergy_path
out_path = args$out_path

# Parameters for determining which biomarker(s) to calculate
lfc_biomarker = as.logical(toupper(args$lfc_biomarker))
auc_biomarker = as.logical(toupper(args$auc_biomarker))
synergy_biomarker = as.logical(toupper(args$synergy_biomarker))
univariate_biomarker = as.logical(toupper(args$univariate_biomarker))
multivariate_biomarker = as.logical(toupper(args$multivariate_biomarker))

# Print some arguments ----
message("INFO: univariate_biomarker is ", univariate_biomarker)
message("INFO: multivariate_biomarker is ", multivariate_biomarker)
message("INFO: lfc_biomarker is ", lfc_biomarker)
message("INFO: auc_biomarker is ", auc_biomarker)
message("INFO: synergy_biomarker is ", synergy_biomarker)

# Create the treatment columns ----
# Need to drop "cell_set"
trt_cols = sig_cols[!grepl("cell_set", sig_cols)]
message("INFO: Trt cols for are  ", paste(trt_cols, collapse = ", "))

# Check if the output directory exists, if not create it
if (!dir.exists(out_path)) {
  dir.create(out_path)
}

# Call the biomarker functions ----
if (univariate_biomarker) {
  if (lfc_biomarker) {
    message("INFO: Creating univariate biomarkers for ", lfc_column)
    create_univariate_biomarker_table(
      in_path = lfc_path,
      out_path = out_path,
      output_file_name = "median_l2fc_univariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = lfc_column,
      depmap_file = bio_file
    )
  }

  if (auc_biomarker) {
    if (file.exists(drc_file)) {
      message("INFO: Creating univariate biomarkers for ", dr_column)

      # Make sure column names are present
      drc = read_data_table(drc_file)
      trt_cols_auc = intersect(trt_cols, colnames(drc))
      message("INFO: Trt cols for auc are  ", paste(trt_cols_auc, sep = ", "))

      create_univariate_biomarker_table(
        in_path = drc_file,
        out_path = out_path,
        output_file_name = "log2_auc_univariate_biomarkers.csv",
        treatment_columns = trt_cols_auc,
        response_column = dr_column,
        depmap_file = bio_file
      )
    } else {
      warning("DRC file does not exist. Skipping univar with AUCs.")
    }
  }

  if (synergy_biomarker) {
    if (file.exists(synergy_path)) {
      message("INFO: Creating univariate biomarkers for ", synergy_col)
      create_univariate_biomarker_table(
        in_path = synergy_path,
        out_path = out_path,
        output_file_name = "synergy_univariate_biomarkers.csv",
        treatment_columns = trt_cols,
        response_column = synergy_col,
        depmap_file = bio_file
      )
    } else {
      warning(paste0("Synergy file not found at ", synergy_file,
                     "\nSkipping univariate biomarkers with synergy scores."))
    }
  }
}

if (multivariate_biomarker) {
  if (lfc_biomarker) {
    message("INFO: Creating multivariate biomarkers for ", lfc_column)
    create_multivariate_biomarker_table(
      in_path = lfc_path,
      out_path = out_path,
      output_file_name = "median_l2fc_multivariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = lfc_column,
      depmap_file = bio_file
    )
  }

  if (auc_biomarker) {
    if (file.exists(drc_file)) {
      message("INFO: Creating multivariate biomarkers for ", dr_column)

      # Make sure column names are present
      drc = read_data_table(drc_file)
      trt_cols_auc = intersect(trt_cols, colnames(drc))
      message("INFO: Trt cols for auc are  ", paste(trt_cols_auc, sep = ", "))

      create_multivariate_biomarker_table(
        in_path = drc_file,
        out_path = out_path,
        output_file_name = "log2_auc_multivariate_biomarkers.csv",
        treatment_columns = trt_cols_auc,
        response_column = dr_column,
        depmap_file = bio_file
      )
    } else {
      warning("DRC file does not exist. Skipping multivar with AUCs.")
    }
  }

  if (synergy_biomarker) {
    if (file.exists(synergy_path)) {
      message("INFO: Creating multivariate biomarkers for ", synergy_col)
      create_multivariate_biomarker_table(
        in_path = synergy_path,
        out_path = out_path,
        output_file_name = "synergy_multivariate_biomarkers.csv",
        treatment_columns = trt_cols,
        response_column = synergy_col,
        depmap_file = bio_file
      )
    } else {
      warning(paste0("Synergy file not found at ", synergy_file,
                     "\nSkipping multivariate biomarkers with synergy scores."))
    }
  }
}
