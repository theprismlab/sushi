options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(magrittr))
source("collapse_replicates/collapse_replicates_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action= "store_true", default= TRUE,
                    help= "Print extra output [default]")
parser$add_argument("-q", "--quietly", action= "store_false", dest= "verbose", 
                    help= "Print little output")
parser$add_argument("-c", "--lfc", default= "l2fc.csv",
                    help= "path to file containing l2fc values")
parser$add_argument("--sig_cols", default= "cell_set,pert_name,pert_dose,pert_dose_unit,day,pert_vehicle",
                    help= "columns used to identify a unique condition")
parser$add_argument("--cell_line_cols", default= "pool_id,depmap_id,lua",
                    help= "Columns that can describe a cell line")
parser$add_argument("--collapsed_l2fc_file", default = "collapsed_l2fc.csv",
                    help = "Name of the file to be stored in the output directory.")
parser$add_argument("--mt_filter", type = "logical", default = TRUE, help = "Filter out outlier treatment pools")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# If the output file already exists, remove it ----
delete_existing_files(args$out, "collapsed_l2fc")

# Collapse biological replicates ----
lfc_values= read_data_table(args$lfc)
sig_cols= unlist(strsplit(args$sig_cols, ","))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))

print("Collapsing biological replicates ...")
collapsed_l2fc= collapse_bio_reps(l2fc= lfc_values, sig_cols= sig_cols, cell_line_cols= cell_line_cols)

# Write out initial collapsed_l2fc file
collapsed_l2fc_outpath = file.path(args$out, "collapsed_l2fc_original.csv")
message("Writing out initial collapsed l2fc file to ", collapsed_l2fc_outpath)
write_out_table(collapsed_l2fc, collapsed_l2fc_outpath)

if (args$mt_filter == TRUE) {
  message("Monotonicity QC filter: On")
  # Set up lists of column names for inputs
  trt_cl_cols = unique(c(setdiff(sig_cols, c("pert_dose", "pert_dose_unit", "pert2_dose", "pert2_dose_unit")),
                         cell_line_cols))
  # trt_cl_cols are used to group all doses of a perturbation for a cell line together
  trt_pool_cols = unique(c(sig_cols, intersect(cell_line_cols, c("cell_set", "pool_id"))))
  # trt_pool_cols are used to group all lines of a pool in each perturbation together

  # trt_cell_set_cols are used to group all lines of a cell set in each perturbation together
  # "cell_set" may already exist in "sig_cols", so only add "cell_set" if it is missing
  if ("cell_set" %in% sig_cols) {
    trt_cell_set_cols = sig_cols
  } else {
    messsage("Attempting to add cell_set to sig_cols ...")
    if ("cell_set" %in% names(collapsed_l2fc)) {
      trt_cell_set_cols = c(sig_cols, "cell_set")
    } else {
      # Error out if "cell_set is not present"
      stop("Cannot find cell_set column in collapsed_l2fc - check sig_cols and/or cell_line_cols!")
    }
  }

  # Call monotonicity check and flag treatment pools
  monotonicity = get_monotonicity(collapsed_l2fc, trt_cl_cols)
  flagged_trt_pools = flag_breaks(monotonicity, trt_pool_cols, trt_cell_set_cols, pool_cutoff = 0.25)

  # Check if any pert + cell sets were flagged
  failed_cell_sets = flagged_trt_pools |>
    dplyr::filter(cell_set_pass_ratio < 0.5) |>
    dplyr::distinct(dplyr::across(tidyselect::all_of(trt_cell_set_cols)))

  # Print out flagged pert + cell sets
  if (nrow(failed_cell_sets) > 0) {
    message("The following ", nrow(failed_cell_sets),
            " cell_sets were dropped due to a large number of pools breaking monotonicity.")
    print(failed_cell_sets)
  }

  # Create a qc output directory if one does not exist
  if (!dir.exists(file.path(args$out, "qc_tables"))) {
    dir.create(file.path(args$out, "qc_tables"))
  }

  # Write out treatment pool qc table as a csv
  outlier_pool_outpath = file.path(args$out, "qc_tables", "trt_pools_qc.csv")
  message("Writing out monotonicity filter file to ", outlier_pool_outpath)
  write_out_table(flagged_trt_pools, outlier_pool_outpath)

  # Filter out failing pools from collapsed_l2fc
  failed_pools = flagged_trt_pools |> dplyr::filter(!is.na(cell_set_flag))
  message("Filtering out ", nrow(failed_pools), " pools from collapsed_l2fc")
  collapsed_l2fc = collapsed_l2fc |> dplyr::anti_join(failed_pools, by = trt_pool_cols)
} else {
  message("Monotonicity QC filter: Off")
}

# Write out file ----
collapsed_l2fc_outpath = file.path(args$out, args$collapsed_l2fc_file)
message("Writing out collapsed l2fc file to ", collapsed_l2fc_outpath)
write_out_table(collapsed_l2fc, collapsed_l2fc_outpath)

# Ensure that collapsed file was successfully generated ----
check_file_exists(collapsed_l2fc_outpath)