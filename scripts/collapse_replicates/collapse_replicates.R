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
parser$add_argument("--mt_filter", type = "logical", default = FALSE, help = "Filter out outlier treatment pools")
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

if (args$mt_filter) {
  message("Monotonicity QC filter: On")
  outlier_pool = get_monotonicity(collapsed_l2fc, cell_line_cols)

  # Create a qc output directory if one does not exist
  if (!dir.exists(file.path(args$out, "qc_tables"))) {
    dir.create(file.path(args$out, "qc_tables"))
  }
  # Write out outlier pool table as a csv
  outlier_pool_outpath = file.path(args$out, "outlier_trt_pools.csv")
  message("Writing out monotonicity filter file to ", outlier_pool_outpath)
  write_out_table(outlier_pool, outlier_pool_outpath)

  message("Filtering out ", nrow(outlier_pool), " pools from collapsed_l2fc")
  failed_pools = outlier_pool |> dplyr::filter(outlier == TRUE)
  join_cols = c("pert_plate", "pert_name", "pert_dose", "cell_set", "pool_id")
  collapsed_l2fc = collapsed_l2fc |> anti_join(outlier_pool, by = join_cols)
}

# Write out file ----
collapsed_l2fc_outpath = file.path(args$out, args$collapsed_l2fc_file)
print(paste0('Writing out collapsed l2fc file to ', collapsed_l2fc_outpath))
write_out_table(collapsed_l2fc, collapsed_l2fc_outpath)

# Ensure that collapsed file was successfully generated ----
check_file_exists(collapsed_l2fc_outpath)