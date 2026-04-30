options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
source("compute_l2fc/compute_l2fc_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                    help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose",
                    help = "Print little output")
parser$add_argument("-c", "--normalized_counts", default = "normalized_counts.csv",
                    help = "path to file containing normalized counts")
parser$add_argument("-ct", "--control_type", default = "ctl_vehicle",
                    help = "pert_type to use as control")
parser$add_argument("--cell_line_cols", default = "cell_set,pool_id,depmap_id,lua",
                    help = "Columns that describe a unique cell line")
parser$add_argument("--sig_cols", default = "pert_name,pert_dose,pert_dose_unit,day",
                    help = "Columns that identify unique perturbations")
parser$add_argument("--ctrl_cols", default = "cell_set,day",
                    help = "Columns used to identify unique control groups")
parser$add_argument("--bio_rep_col", default = "", help = "Column describing the biological replicate.")
parser$add_argument("-ccn", "--count_col_name", default = "log2_normalized_n",
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("-ff", "--filter_failed_lines", type = "logical",
                    help = "Filter out failed cell lines from the output file")
parser$add_argument("--pool_qc_path", default = "", help = "Path to pool level QC file")
parser$add_argument("--plate_cell_qc_path", default = "", help = "Path to pool level QC file")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")
parser$add_argument("--output_file", default = "l2fc.csv", help = "Name of file to output l2fc values")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up parameters ----
normalized_counts = read_data_table(args$normalized_counts)
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
sig_cols = unlist(strsplit(args$sig_cols, ","))
bio_rep_col = args$bio_rep_col
control_type = args$control_type
count_col_name = args$count_col_name
pool_qc_path = args$pool_qc_path
plate_cell_qc_path = args$plate_cell_qc_path

# If l2fc files already exist, remove them ----
delete_existing_files(args$out, "^l2fc")

print("Collapsing tech reps and computing log-fold change ...")
l2fc= compute_l2fc(normalized_counts= normalized_counts,
                   control_type= control_type,
                   sig_cols= sig_cols,
                   bio_rep_col = bio_rep_col,
                   ctrl_cols= ctrl_cols,
                   count_col_name= count_col_name,
                   cell_line_cols= cell_line_cols)

# If filter_failed_lines is TRUE, filter out failed cell lines from the output file ----
if (args$filter_failed_lines) {
  # Check if qc files were provided
  if (pool_qc_path == "") {
    append_critical_output("If filter_failed_lines is TRUE, please provide a path to the pool level QC file.",
                           output_path= args$out)
    stop("If filter_failed_lines is TRUE, please provide a path to the pool level QC file.")
  }
  if (plate_cell_qc_path == "") {
    append_critical_output("If filter_failed_lines is TRUE, please provide a path to the cell line level QC file.",
                           output_path= args$out)
    stop("If filter_failed_lines is TRUE, please provide a path to the cell line level QC file.")
  }

  # Write out the unfiltered l2fc file
  message("Writing out unfiltered l2fc file ...")
  l2fc_unfiltered_outpath = file.path(args$out, "l2fc_original.csv")
  write_out_table(l2fc, l2fc_unfiltered_outpath)

  # Filter failing pools
  pool_qc = read_data_table(pool_qc_path)
  join_cols = c("pert_plate", "pcr_plate", "pcr_well", "cell_set", "pool_id")
  failed_pools = pool_qc |> dplyr::filter(!is.na(well_flag)) |> dplyr::select(all_of(join_cols))

  message("Removing ", nrow(failed_pools), " failing pools ...")
  l2fc = l2fc |> dplyr::anti_join(failed_pools, by = join_cols)

  # Filter failing cell lines
  plate_cell_qc = read_data_table(plate_cell_qc_path)
  join_cols = intersect(c(cell_line_cols, "pert_plate"), colnames(plate_cell_qc))
  failed_lines_pert_plate = plate_cell_qc %>% filter(qc_pass_pert_plate == FALSE) %>% select(all_of(join_cols))

  message("Removing ", nrow(failed_lines_pert_plate), " failing cell lines ...")
  l2fc = l2fc %>% anti_join(failed_lines_pert_plate, by = join_cols)
}

# Write out file ----
l2fc_outpath = file.path(args$out, args$output_file)
print(paste0('Writing out l2fc file to ', l2fc_outpath))
write_out_table(l2fc, l2fc_outpath)

# Ensure that l2fc file was successfully generated ----
check_file_exists(l2fc_outpath)
