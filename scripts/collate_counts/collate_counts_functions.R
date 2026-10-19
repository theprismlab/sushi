#' validate_detected_flowcells
#' 
#' This function checks the table of expected flowcells using the table of detected flowcells.
#' Any flowcells that are expected but are not detected are printed and a warning is printed.
#' 
#' @param detected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
#' @param expected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
validate_detected_flowcells= function(detected_flowcells, expected_flowcells) {
  # Use dplyr::anti_join to filter out row in expected_flowcells that appear in detected_flowcells.
  missing_flowcells= expected_flowcells %>% dplyr::anti_join(detected_flowcells, by= c('flowcell_name', 'flowcell_lane'))

  # Print a warning if there are expected flowcells that were not detected!
  if(nrow(missing_flowcells) != 0) {
    print('WARNING: The following flowcells/lanes specified in the sample meta were not detected in the fastq reads.')
    print(missing_flowcells)
    print('Check that the sample meta is correct or that all fastq files are in the correct directory.')
  }
}

#' collate_fastq_reads
#' 
#' This function takes in the fastq reads (uncollapsed_raw_counts) and
#' filters for reads coming from flowcells specified in the sample meta.
#' The function then sums up the reads across specified sequencing index columns and
#' maps the sequencing index columns to the ID columns.
#' 
#' @param uncollapsed_raw_counts Dataframe of reads from all the fastq files with the following columns -
#'                    "flowcell_name", "flowcell_lane", "index_1", "index_2", "forward_read_barcode", and "n".
#'                    The flowcell columns are optional. If they do not exists, flowcell filters will be skipped.
#' @param sample_meta Sample metadata generate for the project which may contain the following columns - 
#'                    "flowcell_names", "flowcell_lanes", "index_1", "index_2". The sample meta MUST contain
#'                    "flowcell_names" and "flowcell_lanes" for filtering.
#' @param sequencing_index_cols Sequencing columns from the sample meta that the counts should be collapsed on.
#'                              These columns should be a subset of the four sequencing related columns in the
#'                              sample meta - "flowcell_names", "flowcell_lanes", "index_1", and "index_2". They 
#'                              should also uniquely identify every PCR well. This parameter defaults onto 
#'                              the following columns: "index_1", "index_2".
#' @param id_cols ID columns from the sample meta that uniquely identify every PCR well. These columns should not 
#'                include any sequencing related columns. This parameter defaults onto "pcr_plate", "pcr_well". This 
#'                parameter can also be a list of the sample conditions columns as long as they uniquely identify every
#'                PCR well. For example "cell_set", "pert_name", "dose", "day", "bio_rep", "tech_rep" can also be used.
#' @param known_barcodes A vector of known PRISM barcodes. If a read does not match a barcode in this list,
#'                       then its sequence is reassigned to "unknown_reads".
#' @param reverse_index2 Index 2 should be reversed if the sequencer uses a reverse complement workflow. 
#'                       Defaults to FALSE.
#' @param barcode_col String name of the column in uncollapsed_raw_counts that contains the sequences.  
#' @returns Returns a dataframe with columns specified by the id_cols along with barcode_col, and "n".
#' @import tidyverse
#' @import data.table
collate_fastq_reads= function(uncollapsed_raw_counts, sample_meta, 
                              sequencing_index_cols= c('index_1', 'index_2'),
                              id_cols= c('pcr_plate', 'pcr_well'),
                              known_barcodes,
                              reverse_index2= FALSE,
                              barcode_col= 'forward_read_barcode',
                              low_abundance_threshold= 20,
                              out_dir) {
  require(tidyverse)
  require(data.table)
  
  # Validation: Check that sequencing_index_cols exist in the sample_meta ----
  # Error out if a sequencing_index_col is not in the sample_meta.
  if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
    stop('The above sequencing_index_cols are NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols in the sample meta are filled out ----
  # Check for rows in sequencing_index_cols that equate to empty - NA, "NA", "", " "
  # Error out if the sequencing_index_cols are not filled out in the sample meta.
  if(!validate_columns_entries(sequencing_index_cols, sample_meta)) {
    stop('One or more sequencing_index_cols in the sample meta is not filled out.')
  }
  
  # Validation: Check that sequencing_index_cols uniquely identify every row of the sample_meta ----
  # Eror out if a sequencing_index_col does not appear in the sample_meta.
  if(!validate_unique_samples(sequencing_index_cols, sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of sequencing index columns.')
    stop('The specified sequencing index columns do NOT uniquely identify every PCR well.')
  }
  
  # Validation: Check that id_cols exist in the sample meta ----
  # Error out if an id_col is not detected in the sample_meta.
  if(!validate_columns_exist(id_cols, sample_meta)) {
    stop('One or more id_cols is NOT present in the sample meta.')
  }
  
  # Validation: Check that id_cols uniquely identify every row of the sample_meta ----
  if(!validate_unique_samples(id_cols, sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of ID columns.')
    stop('The specified ID columns do NOT uniquely identify every PCR well.')
  }
  
  # Reverse index 2 if specified ----
  if(reverse_index2) {
    if('index_2' %in% colnames(sample_meta)) {
      print("Reverse-complementing index 2 barcode ...")
      sample_meta$index_2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
    } else {
      stop('Reverse index 2 is set to TRUE, but index_2 does not exists.')
    }
  }
  
  # Create sequence map ----
  # Sequencing map is used to map combinations of the sequencing_index_cols to combinations of the id_cols.
  sequencing_map= sample_meta %>% dplyr::distinct(pick(all_of(c(sequencing_index_cols, id_cols))))
  
  # Validation: Check that mapping is one to one ----
  # Make sure that the mapping from sequencing_index_cols to id_cols is 1 to 1.
  # Code below groups on the sequencing_index_cols and filter for combinations that map to more than one id_col combination.
  check_mapping= sequencing_map %>% dplyr::group_by(pick(all_of(sequencing_index_cols))) %>%
    dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
  if(nrow(check_mapping) > 0) {
    print('The following sequencing locations map to multiple conditions.')
    print(check_mapping)
    stop('The sequencing index columns do not map 1 to 1 to the ID columns.')
  }
  
  # Determine if 'flowcell_name' and flowcell_lane' are present in the uncollapsed raw counts file ----
  # If the columns are present, assume that uncollapsed_raw_counts is from Nori and filter for only valid flowcells
  # If not, the file could be from a MiSeq run or something else outside of Nori, so skip this filter step.
  if(base::all(c('flowcell_name', 'flowcell_lane') %in% colnames(uncollapsed_raw_counts))) {
    print('Detecting flowcell_name and flowcell_lane. Filtering for valid flowcells.')
    # Validation: Check that flowcell_names and flowcell_lanes exist in the sample meta ----
    if(!validate_columns_exist(c('flowcell_names', 'flowcell_lanes'), sample_meta)) {
      stop('The above column(s) are NOT present in the sample meta.')
    }
    
    # Determine which flowcell names + lanes are expected ----
    # "flowcell_names" and "flowcell_lanes" are strings that can contain more than one item.
    # Columns can be parsed by splitting on the characters , ; :
    # If there are multiple lane names and lane numbers, this will take the Cartesian product!
    # Note: fread and read.csv keeps commas, read_csv DROPS commas
    expected_flowcells= sample_meta %>% dplyr::distinct(flowcell_names, flowcell_lanes) %>%
      dplyr::mutate(flowcell_name= base::strsplit(flowcell_names, split='[,;:]', fixed=F),
                    flowcell_lane= base::strsplit(flowcell_lanes, split='[,;:]', fixed=F)) %>% 
      tidyr::unnest(cols= flowcell_name) %>% tidyr::unnest(cols= flowcell_lane) %>%
      dplyr::mutate(flowcell_lane= as.numeric(flowcell_lane))
    print(paste0('Expected ', nrow(expected_flowcells), ' unique flowcell + lane combos in the sample meta:'))
    print(expected_flowcells)
    # Note: This code uses base::strsplit and tidyr::unnest from an older version of tidyverse.
    # If there is any update to the tidyverse verision, this can be refactored to use
    # tidyr::separate_longer_delim
    
    # Validation: Check if all expected flowcell name + lanes are detected ----
    # Check that all expected flowcell name + lanes are present in uncollapsed_raw_counts.
    # Print warning if a flowcell is expected but not detected.
    detected_flowcells= uncollapsed_raw_counts %>% dplyr::distinct(flowcell_name, flowcell_lane)
    print(paste0('Identified ', nrow(detected_flowcells), ' unique flowcell + lane combos in the uncollapsed raw counts.'))
    print(detected_flowcells)
    validate_detected_flowcells(detected_flowcells, expected_flowcells)
    
    # Filter for expected flowcells ----
    # Filter using inner join with  merge from data.table instead of dplyr join to improve performance.
    # Dedupe to the join keys first. One flowcell_name + flowcell_lane pair can appear under several
    # flowcell_names/flowcell_lanes groups - a lane 7 listed in both "1,4,7" and "2,5,7,8" - so joining
    # on a table that still carries those columns duplicates every read in a shared lane.
    valid_flowcells= expected_flowcells %>% dplyr::distinct(flowcell_name, flowcell_lane)
    uncollapsed_raw_counts= data.table::merge.data.table(uncollapsed_raw_counts, data.table::setDT(valid_flowcells),
                                                         by= c('flowcell_name', 'flowcell_lane'), allow.cartesian= FALSE)

    # Reads that survive the flowcell filter are the denominator for index purity below. It is recorded
    # here because attaching the flowcell group columns can legitimately expand the read table.
    valid_read_total= sum(uncollapsed_raw_counts$n)

    # Attach the flowcell group columns if they are used as sequencing index columns ----
    # This expansion is undone by the sequencing_map join below, which drops the copies carrying the
    # wrong group. Without flowcell_lanes to tell the groups apart both copies would survive that join,
    # so error out on that combination rather than silently double counting reads.
    carried_flowcell_cols= base::intersect(c('flowcell_names', 'flowcell_lanes'), sequencing_index_cols)
    if(length(carried_flowcell_cols) > 0) {
      flowcell_groups= expected_flowcells %>%
        dplyr::distinct(pick(all_of(c('flowcell_name', 'flowcell_lane', carried_flowcell_cols))))
      shared_lanes= flowcell_groups %>% dplyr::group_by(flowcell_name, flowcell_lane) %>%
        dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
      if(nrow(shared_lanes) > 0 && !'flowcell_lanes' %in% sequencing_index_cols) {
        print('The following lanes are shared by multiple flowcell groups in the sample meta:')
        print(shared_lanes)
        stop('Shared lanes require "flowcell_lanes" in sequencing_index_cols to avoid double counting reads.')
      }
      uncollapsed_raw_counts= data.table::merge.data.table(uncollapsed_raw_counts, data.table::setDT(flowcell_groups),
                                                           by= c('flowcell_name', 'flowcell_lane'),
                                                           allow.cartesian= nrow(shared_lanes) > 0)
    }

  } else {
    print('Flowcell_name and/or flowcell_lane are not detected in raw_counts_uncollapsed.')
    print('Proceeding without filtering flowcells.')
    valid_read_total= sum(uncollapsed_raw_counts$n)
  }

  # Validation: Check that sequencing_index_cols exist in uncollapsed_raw_counts ----
  if(!validate_columns_exist(sequencing_index_cols, uncollapsed_raw_counts)) {
    stop('Some sequencing_index_cols are NOT present in the uncollapsed_raw_counts')
  }
  
  # Create summed_reads file ----
  print('Summing up reads.')
  # Performing inner join with data.table instead of dplyr
  summed_reads= data.table::merge.data.table(uncollapsed_raw_counts, sequencing_map, by= sequencing_index_cols)
  # Code below is checking if a barcode is in the list of known barcodes.
  # If the barcode is not in the list of known barcodes and its counts is below the low_abundance_threshold, 
  # then the barcode is replaced with the string "unknown_low_abundance_barcode".
  # Function := performs the mutate in place without copying the dataframe.
  # Functions fifelse and %chin% are just faster data.table versions of ifelse and %in%.
  summed_reads[, c(barcode_col) := data.table::fifelse(get(barcode_col) %chin% unique(known_barcodes) |
                                                         n >= low_abundance_threshold,
                                                       get(barcode_col), 'unknown_low_abundance_barcode')]

  # Use data.table to group by id_cols and barcode_col and sum up reads across flowcells.
  summed_reads= summed_reads[, .(n= sum(n)), by= c(id_cols, barcode_col)]
  
  # Calculate index purity ----
  # This is only accurate if the Nori input file is small enough to fit into a chunk.
  index_purity= sum(summed_reads$n) / valid_read_total
  # if index_purity is NA, make it 0
  if(is.na(index_purity)){
    index_purity = 0
  }
  
  # Throw an error if the purity is greater than 1.
  # Throw a warning if the purity is below 0.5.
  print(paste0('Index purity in chunk: ', round(index_purity, 4)))
  if(index_purity > 1) {
    stop('ERROR: Chunk index purity is greater than 1!')
  } else if(index_purity < 0.5) {
    print('Warning: Low index purity!')
  } else {}
  
  # Return list of two dfs with known or unknown read counts ----
  print('Completing collate_fastq_reads.')
  return(list(prism_barcode_counts= summed_reads[summed_reads[[barcode_col]] %chin% known_barcodes,], 
              unknown_barcode_counts= summed_reads[!summed_reads[[barcode_col]] %chin% known_barcodes,]))
}

# Record sample_meta rows where either flowcell_names or flowcell_lanes are NA.
check_flowcell_names_in_sample_meta <- function(sample_meta, out_dir) {
  # Check for rows in flowcell_names that equate to empty - NA, "NA", "", " "
  # Error out if the flowcell_names are not filled out in the sample meta.
  missing_flowcell_names= sample_meta %>% dplyr::filter(if_any(all_of(c('flowcell_names')), ~ . %in% c(NA, 'NA', '', ' ')))
  if(nrow(missing_flowcell_names) > 0) {
    append_critical_output('The following rows in the sample meta are not filled out for flowcell_names.', out_dir)
    append_critical_output(missing_flowcell_names, out_dir) # show the empty rows
    append_critical_output('This means that these samples will not have any data in subsequent outputs. Ensure this is intentional.', out_dir)
    append_critical_output('These samples will be saved in the missing_flowcells.csv file', out_dir)
    # Write out missing flowcell names
    write.csv(missing_flowcell_names, paste0(out_dir, '/missing_flowcells.csv'), row.names= FALSE, quote= FALSE)
  } else {
    append_critical_output('All flowcell_names and lanes are filled out in the sample meta.', out_dir)
  }
}