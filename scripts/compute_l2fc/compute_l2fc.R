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
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("-ct", "--control_type", default="ctl_vehicle", help="pert_type to use as control")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--bio_rep_col", default = "", help = "Column describing the biological replicate.")
parser$add_argument("--ctrl_cols", default="cell_set,day",
                    help = "columns used to collapse controls to generate l2fc")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id,lua",
                    help = "Columns that can describe a cell line")
parser$add_argument("-ccn", "--count_col_name", default="log2_normalized_n",
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("--output_file", default = "l2fc.csv", help = "Name of file to output l2fc values")
parser$add_argument("-o","--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-ff", "--filter_failed_lines", type="logical",
                    help = "Filter out failed cell lines from the output file")
parser$add_argument("-qc", "--qc_path", default="", help = "Path to cell line level QC file")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up parameters and run compute_l2fc ----
control_type = args$control_type
normalized_counts= read_data_table(args$normalized_counts)
sig_cols = unlist(strsplit(args$sig_cols, ","))
bio_rep_col = args$bio_rep_col
ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
count_col_name = args$count_col_name
qc_path = args$qc_path

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

# If filter_failed_lines is TRUE, filter out failed cell lines from the output file
if (args$filter_failed_lines) {
  if (qc_path == "") {
    append_critical_output("If filter_failed_lines is TRUE, please provide a path to the cell line level QC file.",
                           output_path= args$out)
    stop("If filter_failed_lines is TRUE, please provide a path to the cell line level QC file.")
  }
  # Write out the unfiltered l2fc file
  print("Writing out unfiltered l2fc file ...")
  l2fc_unfiltered_outpath= paste(args$out, "l2fc_original.csv", sep= "/")
  write_out_table(l2fc, l2fc_unfiltered_outpath)

  # Read in QC file and filter lines that fail for an entire pert_plate
  qc_data = read_data_table(args$qc_path)
  join_cols = intersect(c(cell_line_cols, "pert_plate"), colnames(qc_data))
  failed_lines_pert_plate = qc_data %>% filter(qc_pass_pert_plate == FALSE) %>% select(all_of(join_cols))

  l2fc = l2fc %>% anti_join(failed_lines_pert_plate, by= join_cols)
  # Filter lines that fail for a pcr_plate
  #join_cols = c(cell_line_cols, "pcr_plate")
  #failed_lines_pcr_plate = qc_data %>% filter(qc_pass == FALSE) %>% select(all_of(join_cols))
  #l2fc = l2fc %>% anti_join(failed_lines_pcr_plate, by= join_cols)
}

# Write out file ----
l2fc_outpath = file.path(args$out, args$output_file)
print(paste0('Writing out l2fc file to ', l2fc_outpath))
write_out_table(l2fc, l2fc_outpath)

# Ensure that l2fc file was successfully generated ----
check_file_exists(l2fc_outpath)
