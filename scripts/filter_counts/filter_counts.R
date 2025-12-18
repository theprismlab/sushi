options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(tidyverse)) # load last - after dplyr
suppressPackageStartupMessages(library(httr)) # fread
suppressPackageStartupMessages(library(jsonlite)) # fread
suppressPackageStartupMessages(library(data.table)) # fread
source("filter_counts/filter_counts_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument('--prism_barcode_counts', default= 'prism_barcode_counts.csv', help= 'Path to prism_barcode_counts.csv')
parser$add_argument('--sample_meta', default= 'sample_meta.csv', help= 'Path to sample_meta.csv')
parser$add_argument('--cell_set_and_pool_meta', default= 'cell_set_and_pool_meta.csv', 
                    help= 'Path to cell_set_and_pool_meta.csv')
parser$add_argument('--cell_line_meta', default= 'cell_line_meta.csv', help= 'Path to cell_line_meta.csv')
parser$add_argument('--CB_meta', default= 'CB_meta.csv', help= 'Path to CB_meta.csv')
parser$add_argument('--id_cols', default= 'pcr_plate,pcr_well', 
                    help= 'List of sample_meta column names used to identify every PCR well')
parser$add_argument("--barcode_col", default= "forward_read_barcode",
                    help= "Name of the column in uncollapsed_raw_counts that contains the barcode sequences.")
parser$add_argument("--rm_data", type="logical", help = "Remove bad experimental data")
parser$add_argument("--annot_file", default = "annotated_counts.csv", help = "File for annotated_counts output.")
parser$add_argument("--filt_file", default = "filtered_counts.csv", help = "File for filtered_counts output.")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--screen", default="", help = "Screen")
parser$add_argument("--filter_skipped_wells", type="logical", help = "Filter out skipped wells")
parser$add_argument("--api_key", default="", help = "API key for the clue api")
parser$add_argument("--skipped_wells", default="skipped_wells.csv", help = "Name of the skipped_wells file")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

# If output files already exist, remove them ----
delete_existing_files(args$out, "filtered_counts")
delete_existing_files(args$out, "annotated_counts")

# Read in all input files ----
prism_barcode_counts= read_data_table(args$prism_barcode_counts)
sample_meta= read_data_table(args$sample_meta)
cell_set_and_pool_meta= read_data_table(args$cell_set_and_pool_meta)
cell_line_meta= read_data_table(args$cell_line_meta)
CB_meta= read_data_table(args$CB_meta)
screen= args$screen
api_key= args$api_key

# If the skipped wells file exists, read it in ----
if(file.exists(args$skipped_wells)){
  skipped_wells <- read_data_table(args$skipped_wells)
  skipped_wells <- skipped_wells %>%
    drop_na(pool_id) %>%
    dplyr::rename(
      replicate_plate = replicate,
      pcr_well = assay_well_position
    )
} else {
  print("No skipped wells file found, not filtering skipped wells.")
}

# Convert input strings into vectors ----
id_cols= unlist(strsplit(args$id_cols, ","))

# Remove any duplicate depmap_id + lua + cell_set in cell_set_meta ----
duplicate_ids <- cell_set_and_pool_meta %>%
  dplyr::count(depmap_id, lua, cell_set, name = "count") %>%
  dplyr::filter(count > 1)

if (nrow(duplicate_ids) > 0) {
  print("The following lua + depmap_id combinations are duplicated in a single cell set")
  print(duplicate_ids)
  print("FILTER_COUNTS will continue without considering these duplicated IDs.")
  cell_set_and_pool_meta %<>% dplyr::anti_join(duplicate_ids, by = c("depmap_id", "lua", "cell_set"))
}

# Run filter_raw_reads -----
print('Calling filter_raw_reads ...')
module_outputs= filter_raw_reads(prism_barcode_counts= prism_barcode_counts, 
                                 sample_meta= sample_meta,
                                 cell_set_and_pool_meta= cell_set_and_pool_meta,
                                 cell_line_meta= cell_line_meta,
                                 CB_meta= CB_meta,
                                 id_cols= id_cols,
                                 barcode_col= args$barcode_col)

# Validation: Basic file size check ----
if(sum(module_outputs$filtered_counts$n) == 0) {
  stop('All entries in filtered counts are missing!')
}

# Also check that cell line counts are not all zeros.
filtered_counts= module_outputs$filtered_counts
cl_entries= filtered_counts %>% filter_control_barcodes()
if(sum(cl_entries$n) == 0) {
  stop('All cell line counts are zero!')
}

# If we plan to remove any data, first write a copy of the original unfiltered filtered_counts ----
if(args$rm_data | args$filter_skipped_wells){
  filtered_counts_original <- filtered_counts
  write_out_table(filtered_counts_original, file.path(args$out, "filtered_counts_original.csv"))
}

# Remove data ----
print(paste("rm_data:", args$rm_data))
# Remove data if needed
if(args$rm_data){
  print('rm_data is TRUE, removing supplied data.')
  data_to_remove = read_data_table(file.path(args$out, "data_to_remove.csv"))
  print('Data to remove:')
  print(head(data_to_remove))
  filt_rm <- remove_data(filtered_counts, data_to_remove)

  # re-point to what filtered_counts should be
  filtered_counts <- filt_rm
}

# Filter skipped wells if needed ----
if(exists("skipped_wells") && args$filter_skipped_wells && nrow(skipped_wells) > 0){
  print('filter_skipped_wells is TRUE and the skipped wells file contains entries, removing skipped wells.')
  print('Skipped wells:')
  print(head(skipped_wells))
  filt_skipped <- filter_skipped_wells(filtered_counts, skipped_wells)

  # re-point to what filtered_counts should be
  filtered_counts <- filt_skipped
  rows_removed = nrow(filtered_counts_original) - nrow(filtered_counts)
  paste("Number of skipped wells removed: ", rows_removed)
}

# Write out files ----
annot_out_file = file.path(args$out, args$annot_file)
print(paste('Writing annotated counts to: ', annot_out_file))
write_out_table(module_outputs$annotated_counts, annot_out_file)

filtrc_out_file = file.path(args$out, args$filt_file)
print(paste("Writing filtered counts csv to: ", filtrc_out_file))
write_out_table(filtered_counts, filtrc_out_file)

# Ensure that files were successfully generated ----
check_file_exists(annot_out_file)
check_file_exists(filtrc_out_file)
