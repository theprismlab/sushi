options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
source("collate_counts/collate_counts_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument('--raw_counts_uncollapsed', default= "raw_counts_uncollapsed.csv.gz",
                    help= "path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help= "Cell line metadata")
parser$add_argument("--CB_meta", default= "CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument('--sequencing_index_cols', default= 'flowcell_names,index_1,index_2', 
                    help= 'List of sequencing columns in the sample meta.')
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well", 
                    help = "Columns that identify a unique PCR well")
parser$add_argument("--reverse_index2", type="logical", default=FALSE,
                    help= "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("--barcode_col", default= "forward_read_barcode",
                    help= "Name of the column in uncollapsed_raw_counts that contains the barcode sequences.")
parser$add_argument('--low_abundance_threshold', default= 20, 
                    help= 'For unknown barcodes, counts below this threshold will be marked as an unknown barcode.')
parser$add_argument('--chunk_size', default= 10000000, help= 'Integer number of rows for a chunk.')
parser$add_argument("--prism_barcode_counts_file", default = "prism_barcode_counts.csv",
                    help = "Name of file for prism barcodes.")
parser$add_argument("--unknown_barcode_counts_file", default = "unknown_barcode_counts.csv",
                    help = "Name of file for unknown barcodes.")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# If the output files already exist, remove them ----
delete_existing_files(args$out, "barcode_counts")

# Read in metadata files as data.table objects ----
cell_line_meta= read_data_table(args$cell_line_meta)
CB_meta= read_data_table(args$CB_meta)
sample_meta= read_data_table(args$sample_meta)
## if flowcell_lanes one lane, the integer gets read in as numeric, so force this to be a string
if(validate_columns_exist("flowcell_lanes", sample_meta)){
  sample_meta %<>% mutate(flowcell_lanes = as.character(flowcell_lanes))
}

# Check that sample_meta has complete flowcell names and lanes
check_flowcell_names_in_sample_meta(sample_meta, args$out)

# Parse some parameters into vectors ----
sequencing_index_cols= unlist(strsplit(args$sequencing_index_cols, ","))
id_cols= unlist(strsplit(args$id_cols, ","))

# Validation: Check that sequencing_index_cols present in the sample meta ----
if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
  stop('One or more sequencing_index_cols is NOT present in the sample meta.')
}

# Validation: Check that id_cols are present in the sample meta ----
if(!validate_columns_exist(id_cols, sample_meta)) {
  stop('One or more id_cols is NOT present in the sample meta.')
}

# Validation: Check that barcode_col is present in the CB_meta ----
if(!args$barcode_col %in% colnames(CB_meta)) {
  stop('barcode_col is NOT present in the CB_meta.')
}

# Validation: Check that barcode_col is present in the cell_line_meta ----
if(!args$barcode_col %in% colnames(cell_line_meta)) {
  stop('barcode_col is NOT present in the cell_line_meta.')
}

# Validation: Check that barcode_col is present in the raw_counts_uncollapsed ----
raw_counts_uncollapsed_header <- read_data_table(args$raw_counts_uncollapsed, nrows = 0)
if (!args$barcode_col %in% colnames(raw_counts_uncollapsed_header)) {
  stop('barcode_col is NOT present in the raw_counts_uncollapsed.')
}

# Run collate_fastq_reads on chunks of raw_counts_uncollapsed.csv ----
# raw_counts_uncollapsed could be too large to read into memory,
# so collate_fastq_reads is performed on chunks of the raw_counts_uncollapsed file.
chunked_results= process_in_chunks(large_file_path= args$raw_counts_uncollapsed, 
                                   chunk_size= base::strtoi(args$chunk_size), 
                                   action= collate_fastq_reads,
                                   # Parameters for collate_fastq_reads
                                   sample_meta= sample_meta, 
                                   sequencing_index_cols= sequencing_index_cols,
                                   id_cols= id_cols,
                                   known_barcodes= unique(c(cell_line_meta[[args$barcode_col]], 
                                                            CB_meta[[args$barcode_col]])),
                                   reverse_index2= args$reverse_index2,
                                   barcode_col= args$barcode_col,
                                   low_abundance_threshold= as.numeric(args$low_abundance_threshold),
                                   out_dir = args$out)

# From each chunk, extract prism_barcode_counts or unknown_barcode_counts and bind those rows together.
# Then use data.table to aggregate and sum up reads across the chunks.
# data.table functions are faster and less memory intensive.
prism_barcode_counts= data.table::rbindlist(lapply(chunked_results, function(x) x$prism_barcode_counts))
prism_barcode_counts= prism_barcode_counts[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

unknown_barcode_counts= data.table::rbindlist(lapply(chunked_results, function(x) x$unknown_barcode_counts))
unknown_barcode_counts= unknown_barcode_counts[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

# Validation: Basic file size check ----
if(nrow(prism_barcode_counts) == 0) {
  stop('ERROR: Empty file generated. No rows in prism_barcode_counts output.')
}

# Write out files ----
out_file = file.path(args$out, args$prism_barcode_counts_file)
print(paste("Writing prism_barcode_counts.csv to ", out_file))
write_out_table(prism_barcode_counts, out_file)

# Ensure that files were successfully generated ----
check_file_exists(out_file)

out_file = file.path(args$out, args$unknown_barcode_counts_file)
print(paste("Writing unknown_barcode_counts.csv to ", out_file))
write_out_table(unknown_barcode_counts, out_file)

# Ensure that files were successfully generated ----
check_file_exists(out_file)
