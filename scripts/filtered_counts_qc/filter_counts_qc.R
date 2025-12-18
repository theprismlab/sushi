options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales)) # for out of bound handling in plots
suppressPackageStartupMessages(library(ggpmisc)) # with ggplot to add linear fit labels
suppressPackageStartupMessages(library(WGCNA)) # for faster correlations
source("filter_counts_qc/filter_counts_qc_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()

# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument('--raw_counts_uncollapsed', default= "raw_counts_uncollapsed.csv",
                    help= 'Path to file containing uncollapsed raw counts')
parser$add_argument('--chunk_size', default= 10000000, help= 'Integer number of rows for a chunk.')
parser$add_argument('--prism_barcode_counts', default= "prism_barcode_counts.csv", help= 'Path to prism_barcode_counts.csv')
parser$add_argument('--unknown_barcode_counts', default= "unknown_barcode_counts.csv", 
                    help= 'Path to unknown_barcode_counts.csv')
parser$add_argument('--annotated_counts', default= 'annotated_counts.csv', help= 'Path to annotated_counts.csv')
parser$add_argument('--normalized_counts', default= 'normalized_counts.csv', help= 'Path to normalized_counts.csv')
parser$add_argument('--lfc', default= 'l2fc.csv', help= 'Path to l2fc.csv')
parser$add_argument('-s', '--sample_meta', default= 'sample_meta.csv', help= 'Path to sample_meta.csv')
parser$add_argument("--barcode_col", default= "forward_read_barcode",
                    help= "Name of the column in uncollapsed_raw_counts that contains the barcode sequences.")
parser$add_argument('--id_cols', default= 'pcr_plate,pcr_well', help= 'Sample meta columns used to identify every PCR well')
parser$add_argument('--cell_line_cols', default= 'DepMap_ID', help= 'Sushi columns used to identify a read')
parser$add_argument('--sig_cols', default= 'cell_set,pert_name,pert_dose,pert_dose_unit,day',
                    help= 'Sample meta columns used to identify unique treatment conditions')
parser$add_argument('--control_type', default= 'negcon', help= 'Value used in pert_type column to denote negative controls')
parser$add_argument('--count_threshold', default= 40, help= 'Low counts theshold used in some plots')
parser$add_argument('--reverse_index2', type= "logical", default= FALSE,
                    help= 'Switch to reverse complement index_2 for some sequencers')
parser$add_argument('-o', '--out', default= '', help= 'Output path, defaults to working directory')
parser$add_argument("--pseudocount", default=20, help = "pseudo raw count used in normalization")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}

# Read in input files ----
sample_meta= read_data_table(args$sample_meta)
prism_barcode_counts= read_data_table(args$prism_barcode_counts)
unknown_barcode_counts= read_data_table(args$unknown_barcode_counts)
annotated_counts= read_data_table(args$annotated_counts)
if(file.exists(args$normalized_counts)) {
  print('Detecting normalized_counts.')
  normalized_counts= read_data_table(args$normalized_counts)
} else {
  print('NOT finding normalized_counts ...')
  normalized_counts= NA
}
l2fc= read_data_table(args$lfc)

# Parse some input parameters ----
id_cols= unlist(strsplit(args$id_cols, ','))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ','))
sig_cols= unlist(strsplit(args$sig_cols, ','))

# Call QC images function ----
print("Calling QC images ...")
QC_images(raw_counts_uncollapsed_path= args$raw_counts_uncollapsed,
          chunk_size= base::strtoi(args$chunk_size),
          prism_barcode_counts= prism_barcode_counts, 
          unknown_barcode_counts= unknown_barcode_counts,
          annotated_counts= annotated_counts, 
          normalized_counts= normalized_counts, 
          l2fc= l2fc, 
          sample_meta= sample_meta,
          id_cols= id_cols, 
          barcode_col= args$barcode_col,
          cell_line_cols= cell_line_cols,
          sig_cols= sig_cols,
          control_type= args$control_type, 
          count_threshold= as.numeric(args$count_threshold), 
          reverse_index2= args$reverse_index2,
          out= args$out,
          pseudocount = as.numeric(args$pseudocount))
