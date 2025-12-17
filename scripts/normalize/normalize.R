options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("normalize/normalize_functions.R")
source("utils/kitchen_utensils.R")
source("qc_tables/qc_tables_functions.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--norm_method", default = "normalize_then_shift", help = "Normalization method to use")
parser$add_argument("--id_cols", default="pcr_plate,pcr_well",
                    help = "Columns to identify each PCR well")
parser$add_argument("--CB_meta", default="CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument("--read_detection_limit", default = "10",
                    help = "Smallest read count value. Used to compute pseudovalue.")
parser$add_argument("--negcon_cols", default = "pcr_plate,pert_vehicle",
                    help = "List of columns in filtered counts that describe a negative control condition.")
parser$add_argument("--negcon_type", default = "ctl_vehicle",
                    help = "String in the column pert_type that indicates a negative control.")
parser$add_argument("--pseudocount", default = 0, help = "Pseudocount used in normalization.")
parser$add_argument("--req_negcon_reps", default = 6, help = "Number of required negcons per plate.")
parser$add_argument("--output_file", default = "normalized_counts.csv", help = "File name for normalized counts.")
parser$add_argument("-o", "--out", default=getwd(), help= "Output path. Defaults to working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up inputs
norm_method = args$norm_method
filtered_counts = read_data_table(args$filtered_counts)
CB_meta = read_data_table(args$CB_meta)
input_pseudocount = as.integer(args$pseudocount)
input_id_cols = unlist(strsplit(args$id_cols, ","))
read_detection_limit = as.integer(args$read_detection_limit)
negcon_cols = unlist(strsplit(args$negcon_cols, split = ","))
negcon_type = args$negcon_type
output_file = args$output_file
req_negcon_reps = as.integer(args$req_negcon_reps) # Using default value

# Remove normalized files that already exist in the directory
delete_existing_files(args$out, "normalized_counts")

# Identify CBs in PCR wells that can be normalized.
message("Identifying control barcodes to normalize ...")
cb_annots = flag_control_bcs(filtered_count = filtered_counts,
                             CB_meta = CB_meta,
                             id_cols = input_id_cols,
                             negcon_type = negcon_type,
                             cb_mad_cutoff = 1,
                             req_negcon_reps = req_negcon_reps)

# Extract just the necessary columns from cb_annots
cb_annots = cb_annots |> dplyr::select(all_of(c(input_id_cols, "cb_ladder", "cb_name", "keep_cb")))

# Error out if an incorrect normalization method is supplied
if (!norm_method %in% c("normalize_with_pseudocount", "normalize_then_shift")) {
  stop("'", norm_method, "' is not a valid normalization method. 
       Please set NORM_METHOD to either 'normalize_with_pseudocount' or 'normalize_then_shift'.")
}

# Normalize with pseudocount
message("Normalization method: ", norm_method)
message("Pseudocount value: ", input_pseudocount)
normalized_counts = normalize(X = filtered_counts, cb_annots = cb_annots,
                              id_cols = input_id_cols, pseudocount = input_pseudocount)

# Check that an appropriate pseudocount is provide if normalizing with a pseudocount
if (norm_method == "normalize_with_pseudocount" && input_pseudocount < read_detection_limit) {
  message(sprintf("The pseudocount value (%d) is too low for this normalization method.", input_pseudocount))
  warning("Use a pseudocount value greater than the read detection limit of ", read_detection_limit,
          " or use the other normalization method.")
}

# If normalization_then_offset is selected, add offset to normalized values
if (norm_method == "normalize_then_shift") {
  # Verify that there are enough negative controls
  negcon_reps = filtered_counts[pert_type == negcon_type] |>
    dplyr::distinct(across(all_of(input_id_cols))) |>
    dplyr::group_by(pcr_plate) |>
    dplyr::summarise(num_negcon_reps = dplyr::n(), .groups = "drop")

  if (min(negcon_reps$num_negcon_reps) >= req_negcon_reps) {
    # Run pseudovalue addition if enough negcon reps are present
    message("Adding pseudovalue corresponding to a read count of ", read_detection_limit, " ...")
    normalized_counts = add_pseudovalue(normalized_counts, negcon_cols = negcon_cols,
                                        read_detection_limit = read_detection_limit,
                                        negcon_type = negcon_type)
  } else {
    # Error out if no pseudocount is provided and there are not enough negative control replicates
    message("There are not enough negative control replicates on every PCR plate to calculate the pseudovalue.")
    stop("Not enough negative control replicates. Try normalize_with_pseudocount instead.")
  }
}

# Write out file
normcounts_outpath = file.path(args$out, output_file)
message("Writing normalized count file to: ", normcounts_outpath)
write_out_table(normalized_counts, normcounts_outpath)

# Ensure that normalized file was sucessfully generated
check_file_exists(normcounts_outpath)