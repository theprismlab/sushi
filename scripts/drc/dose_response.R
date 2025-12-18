options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("drc/dose_response_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("--l2fc", default="l2fc.csv", help="file containing log fold change values")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id,lua",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--l2fc_column", default="l2fc",
                    help = "column containing log fold change values")
parser$add_argument("--type_col", default="pert_type",
                    help = "column containing perturbation type, only those with trt_cp will be considered.")
parser$add_argument("--screen_type", default = "", help = "Type of PRISM screen.")
parser$add_argument("--cap_for_viability", default = 1.5,
                    help = paste0("The upper threshold for the viability values for before fitting the curves.
                                   Default is 1.5, any viability value above this value will be made equal to 1.5"))
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--out_dir", default= "", help = "Path to the output directory")

# Get command line options, if help option encountered p3rint help and exit
args <- parser$parse_args()

# Read in l2fc file ----
l2fc= read_data_table(args$l2fc)

# Parse some parameters ----
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
sig_cols= unlist(strsplit(args$sig_cols, ","))
l2fc_col= args$l2fc_col
type_col= args$type_col
cap_for_viability= as.numeric(args$cap_for_viability)
build_dir = args$build_dir

# Pull out dose columns from sig_cols
dose_cols = sig_cols[grepl("_dose$", sig_cols)]
print(paste0("Detecting the following dose column(s): ", paste(dose_cols, collapse = ", ")))
if (length(dose_cols) > 2) {
  # Error out if more than two dose columns are detected.
  stop("More than two dose columns were detected!")
}

# create_drc_table is called for each dose column detected.
# For MTS, the function should just be called once.
# For CPS, the function should be called twice.
drc_outputs = list()

for (idx in seq_along(dose_cols)) {
  cur_dose_col = dose_cols[idx]
  print(paste0("Working on DRC over ", cur_dose_col))

  # Filter l2fc dt by dropping NAs in the current dose column
  subset_l2fc = l2fc[!is.na(get(cur_dose_col)), ]

  # Move to the next iteration if subset at current dose column is empty.
  if (nrow(subset_l2fc) == 0) {
    warning(paste0("No avaliable doses detected for ", cur_dose_col, ". Skipping this dose."))
    next()
  }

  # Move to the next iteration if there are very few unique dose values.
  if (length(na.omit(unique(subset_l2fc[[cur_dose_col]]))) <= 4) {
    warning(paste0("Detected four or fewer unique doses in ", cur_dose_col, ". Skipping this dose."))
    next()
  }

  # DRC function should be called once for MTS, EPS, and CPS single agents.
  # DRC function should be called twice for CPS 5x5s.
  # Both pert_dose and pert_dose_unit will be dropped when grouping
  drc_outputs[[idx]] = create_drc_table(
    LFC = subset_l2fc,
    screen_type = args$screen_type,
    cell_line_cols = cell_line_cols,
    treatment_cols = sig_cols[!grepl(cur_dose_col, sig_cols)],
    dose_col = cur_dose_col,
    l2fc_col = l2fc_col,
    cap_for_viability = cap_for_viability
  )
}

dose_response = dplyr::bind_rows(drc_outputs)

# Validation: Check that dose_response is not empty ----
if (nrow(dose_response) == 0) {
  warning("Dose response table is empty. An output file will not be created.")
} else {
  # Write out the DRC table ----
  # Check if the output directory exists, if not create it
  if (!dir.exists(args$out_dir)) {
    dir.create(args$out_dir)
  }

  drc_outpath = file.path(args$out_dir, "DRC_TABLE.csv")
  paste0("Writing DRC_TABLE.csv to ", drc_outpath)
  write_out_table(dose_response, drc_outpath)

  # Check to make sure that the file was generated
  check_file_exists(drc_outpath)
}
