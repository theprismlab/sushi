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
parser$add_argument("-o", "--out", default= getwd(), help= "Output path. Default is working directory")


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

# Write out file ----
collapsed_l2fc_outpath = file.path(args$out, args$collapsed_l2fc_file)
print(paste0('Writing out collapsed l2fc file to ', collapsed_l2fc_outpath))
write_out_table(collapsed_l2fc, collapsed_l2fc_outpath)

# Ensure that collapsed file was successfully generated ----
check_file_exists(collapsed_l2fc_outpath)