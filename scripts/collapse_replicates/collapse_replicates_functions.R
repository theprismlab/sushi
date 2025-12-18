#' validate_num_bio_reps
#' 
#' Function that checks if biological replicates were collapsed by comparing the 
#' number of unique bio_rep annotations (for example 1, 2, 3 or A, B, C) to the
#' maximum number of biological replicates that were collapsed into a value.
#' 
#' @param num_unique_bio_reps A dataframe with the columns "flowcell_name" and "flowcell_lane".
#' @param max_bio_rep_count A dataframe with the columns "flowcell_name" and "flowcell_lane".
validate_num_bio_reps= function(num_unique_bio_reps, max_bio_rep_count) {
  if(num_unique_bio_reps > 1 & max_bio_rep_count == 1) {
    print('Warning - Detecting unique bio_rep annotations, but each cell line + condition only has one biological replicate.')
    print('Check the sample meta and the l2fc file to make sure this is the intended behavior!')
  } else if(num_unique_bio_reps < max_bio_rep_count) {
    stop('Bio_reps were incorrectly collapsed resulting in more replicates than specified.')
  } else if(num_unique_bio_reps > 1 & num_unique_bio_reps > max_bio_rep_count) {
    print('Warning - Number of replicates that were collapses is smaller than the number of unique bio_rep annotations.')
    print('This could be due to a problem in processing or from poorer data/sequencing quality.')
  } else {}
}

#' collapse_bio_reps
#' 
#' Collapses the l2fc values of biological replicates for treatment conditions and
#' computes the MAD/sqrt(n).
#'
#' @param l2fc Dataframe of l2fc values The following columns are required -
#'              depmap_id, lua, counts_flag, mean_n, mean_normalized_n, and l2fc.
#' @param sig_cols List of columns that define an individual condition. This should not include any replicates.
#'                  The columns in this list should be present in the l2fc dataframe.
#' @param cell_line_cols List of columns that define a cell line. Defaults to "project_code" and "depmap_id"
#' @returns - collapsed_counts
collapse_bio_reps= function(l2fc, sig_cols, cell_line_cols= c('depmap_id', 'lua', 'pool_id')) {
  # Validation: Check that sig_cols are present in l2fc ----
  if(validate_columns_exist(sig_cols, l2fc) == FALSE) {
    print(sig_cols)
    stop('Not all sig_cols (printed above) are present in the l2fc file.')
  }
  
  # Validation: Check that cell_line_cols are present in l2fc ----
  if(validate_columns_exist(cell_line_cols, l2fc) == FALSE) {
    print(cell_line_cols)
    stop('Not all cell_line_cols (printed above) are present in the l2fc file.')
  }

  # Collapse bio replicates and calculate median l2fcs
  median_cols = c("l2fc", "l2fc_uncorrected")
  collapsed_counts = l2fc |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(c(cell_line_cols, sig_cols)))) |>
    dplyr::summarise(dplyr::across(tidyselect::any_of(median_cols), median, .names = "median_{.col}"),
                     num_bio_reps = dplyr::n(), .groups = "drop")

  # Validation: Check that replicates were collapsed ----
  if('bio_rep' %in% colnames(l2fc)) {
    num_unique_bio_reps= length(unique(l2fc$bio_rep))
    max_bio_rep_count= max(unique(collapsed_counts$num_bio_reps))
    validate_num_bio_reps(num_unique_bio_reps, max_bio_rep_count)
  }
  
  return(collapsed_counts)
}
