library(argparse)
library(tidyverse)
source("join_metadata/join_metadata_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument('--sample_meta', default= 'sample_meta.csv', help= 'Sample meta data for the sequencing run.')
parser$add_argument('--lfc', default= 'l2fc.csv', help= 'L2FC data.') # level 4
parser$add_argument('--collapsed_lfc', default= 'collapsed_l2fc.csv', help= 'Collapsed l2fc data.') # level 5
parser$add_argument('--sig_cols', default= 'cell_set,pert_name,pert_dose,pert_dose_unit,day',
                    help= 'Columns that uniquely identify a condition.') 
parser$add_argument('--out', default= getwd(), help= 'Path to the output directory.')

args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# Read in files and prepare some parameters ----
sample_meta= read_data_table(args$sample_meta)
sig_cols= unlist(strsplit(args$sig_cols, ","))

# Add sample meta to l2fc table ----
if(file.exists(args$lfc)) {
  print('Attempting to add sample_meta to l2fc file.')
  l2fc= read_data_table(args$lfc)
  
  # Add sample meta columns to l2fc
  if('bio_rep' %in% colnames(sample_meta) & 'bio_rep' %in% colnames(l2fc)) {
    input_cols= c(sig_cols, 'bio_rep')
  } else {
    input_cols= sig_cols
    print('WARNING: No "bio_rep" column detected. Proceeding with just sig_cols.')
  }
  l2fc_with_meta_columns= join_metadata(input_df= l2fc, metadata= sample_meta, key_cols= input_cols)
  
  # Write out
  outpath= file.path(args$out, "l2fc_with_meta_columns.csv")
  print(paste("Writing l2fc_with_meta_columns.csv to ", outpath))
  write_out_table(l2fc_with_meta_columns, outpath)
} else {
  print('WARNING: l2fc.csv does not exist. Skipping this file.')
}

# Add sample meta to collapsed_l2fc table ----
if(file.exists(args$collapsed_lfc)) {
  print('Attempting to add sample_meta to collapsed l2fc.')
  collapsed_l2fc= read_data_table(args$collapsed_lfc)
  
  # Add sample meta columns to collapsed l2fc
  collapsed_l2fc_with_meta_columns= join_metadata(input_df= collapsed_l2fc, metadata= sample_meta, 
                                                  key_cols= sig_cols)
  
  # Write out
  outpath= file.path(args$out, "collapsed_l2fc_with_meta_columns.csv")
  print(paste("Writing collapsed_l2fc_with_meta_columns.csv to ", outpath))
  write_out_table(collapsed_l2fc_with_meta_columns, outpath)

} else {
  print('WARNING: collapsed_l2fc.csv does not exist. Skipping this file.')
}
