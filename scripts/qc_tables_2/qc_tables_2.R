options(cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(argparse)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(jsonlite)
  library(dplyr)
})
source("qc_tables_2/qc_tables_2_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument(
  "--l2fc",
    default = "l2fc.csv", help = "l2fc file"
)
parser$add_argument(
  "-o", "--out",
  default = getwd(), help = "Output path. Default is working directory"
)


args <- parser$parse_args()

# Read in files as data.table objects ----
print(paste0("Reading in ", args$l2fc, "....."))
l2fc <- read_data_table(args$l2fc)

# Check if the output directory exists, if not create it
if (!dir.exists(paste0(args$out, "/qc_tables"))) {
  dir.create(paste0(args$out, "/qc_tables"))
}

# Pool replicate concordance
pool_delta_df <- compute_pool_delta_df(l2fc, 2)

# WRITE OUT RESULTS --------

# Write pool_delta_df  table
pool_delta_outpath <- paste0(args$out, "/qc_tables/pool_delta_table.csv")
print(paste0("Writing out pool_delta_df to ", pool_delta_outpath))
write_out_table(table = pool_delta_df, path = pool_delta_outpath)
check_file_exists(pool_delta_outpath)

paste0("QC module post-lfc completed.")
