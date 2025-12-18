#' compute_l2fc
#' 
#' takes normalized counts and computes log-fold change values as compared to the designated control condition
#'
#' @param normalized_counts - table with normalized_n column and pert_type column that designates the
#'          the control sample
#' @param control_type - string that denotes which samples to compute log fold change against. Matches pert_type field.
#'          negcon by default.
#' @param sig_cols - a vector of column names denoting which values specify each individual signature
#'                    cell_set,pert_name,dose,dose_unit,day by default.
#' @param ctrl_cols - a vector of column names denoting which values specify each individual control condition
#'                    cell_set,day by default.
#' @param count_col_name - a string containing the name of the column to use as counts to calculate l2fc values. 
#'          Generally log2_normalized_n if running on normalized_counts or n if running on filtered_counts
#' @param cell_line_cols - Vector of columns that define a cell line. Defaults to lua, depmap_id, and pool_id
#' @return - l2fc data.frame with l2fc column
#' @export
compute_l2fc= function(normalized_counts,
                       control_type = "ctl_vehicle",
                       sig_cols=c('cell_set','pert_name','pert_dose','pert_dose_unit','day'),
                       bio_rep_col = "bio_rep",
                       ctrl_cols= c('cell_set', 'day'), # will probably be a subset of sig_cols
                       count_col_name="log2_normalized_n",
                       cell_line_cols= c('lua', 'depmap_id', 'pool_id')) {
  
  # Validation: Check that sig_cols are in normalized_counts ----
  if(validate_columns_exist(sig_cols, normalized_counts) == FALSE) {
    print(sig_cols)
    stop('Not all sig_cols (printed above) are present in normalized_counts.')
  }
  
  # Validation: Check that cell_line_cols are in normalized_counts ----
  if(validate_columns_exist(cell_line_cols, normalized_counts) == FALSE) {
    print(cell_line_cols)
    stop('Not all cell_line_cols (printed above) are present in normalized_counts.')
  }
  
  # Collapsing technical replicates ----
  # Detect bio_rep column to be used to collapse technical replicates
  if(bio_rep_col == "") {
    message("No bio_rep column specified. Assuming that there are NO biological replicates.")
    sig_bio_rep_cols = sig_cols
  } else {
    sig_bio_rep_cols = unique(c(sig_cols, bio_rep_col))
  }
  
  # collapse tech reps
  grouping_cols = unique(c(cell_line_cols, "pert_type", sig_bio_rep_cols, ctrl_cols))
  message("Collapsing technical replicates on the following columns: ",
          paste(grouping_cols, collapse = ", "))

  collapsed_tech_rep= normalized_counts %>%
    filter_control_barcodes() %>%
    dplyr::filter(!(pert_type %in% c('empty', '', 'CB_only', NA))) %>%
    dplyr::group_by(across(all_of(grouping_cols))) %>%
    dplyr::summarise(mean_n= mean(n),
                     mean_normalized_n= mean(2^(!!rlang::sym(count_col_name)))) %>% dplyr::ungroup()

  # Print out the occurrence of each count of tech_reps
  # print('Number of technical replicate collapsed across all cell lines and biological replicates:')
  # print(collapsed_tech_rep %>% dplyr::group_by(num_tech_reps) %>%
  #         dplyr::summarise(count= dplyr::n()) %>% dplyr::ungroup())

  # Pull out negative controls and collapse any biological replicates ----
  print('Collapsing control conditions on the following columns: ')
  print(unique(c(cell_line_cols, ctrl_cols)))
  controls= collapsed_tech_rep %>% dplyr::filter(pert_type== control_type) %>% 
    dplyr::group_by(pick(all_of(union(cell_line_cols, ctrl_cols)))) %>%
    dplyr::summarise(control_median_n= median(mean_n),
                     control_median_normalized_n = median(mean_normalized_n),
                     num_ctrl_bio_reps = dplyr::n()) %>% dplyr::ungroup()
  
  # Validation: Check that negative controls were extracted ----
  if(nrow(controls) == 0) {
    stop("No samples found for the indicated control_type.")
  }

  # Join neg_cons and compute l2fc ----
  l2fc= collapsed_tech_rep %>% dplyr::filter(!pert_type %in% c(control_type, 'day_0')) %>% 
    dplyr::inner_join(controls, by= union(cell_line_cols, ctrl_cols), relationship= 'many-to-one') %>%
    dplyr::mutate(l2fc = log2(mean_normalized_n) - log2(control_median_normalized_n)) %>%
    dplyr::select(-mean_n, -control_median_n)

  return(l2fc)
}
