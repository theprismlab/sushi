#' Calculate purity metrics
#' 
#' Create the qc table with index purity and cell line purity.
#' 
#' @param raw_counts_uncollapsed_path Path to the nori raw_counts_uncollapsed file.
#' @param unknown_barcode_counts Dataframe of unknown barcodes.
#' @param prism_barcode_counts Dataframe of prism barcodes extracted from the nori. 
#' @param filtered_counts Filtered counts dataframe created from filter_raw_reads.
#' @param value_col String name of the counts column present all the four input dataframes.
#' @param file_path Location to write out the output.
#' @returns Writes out a QC_table to the file_path.
create_qc_table= function(raw_counts_uncollapsed_path, unknown_barcode_counts, 
                          prism_barcode_counts, filtered_counts,
                          value_col= 'n', chunk_size= 10^6, output_path) {
  # Validation: Check that the file at the path exists
  if(!file.exists(raw_counts_uncollapsed_path)) {
    stop('Cannot find the raw counts uncollapsed file.')
  }
  
  # Pull out only the headers of the large file for validation
  rcu_headers= read_data_table(raw_counts_uncollapsed_path, nrows= 0)
  
  # Validation: Check that value_col exists in raw_counts_uncollapsed
  if(!validate_columns_exist(value_col, rcu_headers)) {
    stop(paste0('The column ', value_col, ' was not detected in uncollapsed raw counts.'))
  }
  
  # Validation: Check that value_col exists in unknown_barocde_counts
  if(!validate_columns_exist(value_col, unknown_barcode_counts)) {
    stop(paste0('The column ', value_col, ' was not detected in unknown_barcode_counts.csv'))
  }
  
  # Validation: Check that value_col exists in prism_barocde_counts
  if(!validate_columns_exist(value_col, prism_barcode_counts)) {
    stop(paste0('The column ', value_col, ' was not detected in prism_barcode_counts.csv'))
  }
  
  # Validation: Check that value_col exists in filtered_counts
  if(!validate_columns_exist(value_col, filtered_counts)) {
    stop(paste0('The column ', value_col, ' was not detected in filtered_counts.csv'))
  }
  
  # Determine total number of reads in raw_counts_uncollapsed using chunking
  chunk_sum= process_in_chunks(large_file_path= raw_counts_uncollapsed_path, 
                               chunk_size= chunk_size, 
                               action= function(x) data.table::as.data.table(sum(x[[value_col]])))
  total_num_reads= sum(unlist(chunk_sum))
  
  # Determine number of reads that mapped to valid PCR locations
  # These reads have the correct index barcodes
  total_valid_pcr_reads= sum(unknown_barcode_counts[[value_col]]) + sum(prism_barcode_counts[[value_col]])
  
  # Calculate purities
  # Index purity is the fraction of reads that mapped to valid PCR locations out of the total number of reads.
  index_purity= total_valid_pcr_reads / total_num_reads
  print(paste0('Index purity: ', round(index_purity, 4)))
  # Cell line purity is the fraction of reads that are identified as cell lines or control barcodes out of valid PCR reads.
  cell_line_purity= sum(filtered_counts[[value_col]]) / total_valid_pcr_reads
  print(paste0('Cell line purity: ', round(cell_line_purity, 4)))
  
  # Write out QC table
  qc_table= data.frame(index_purity= index_purity, cell_line_purity= cell_line_purity)
  print(paste0('Writing QC table out to ', output_path))
  qc_table %>% write.csv(output_path, row.names= FALSE, quote= FALSE)
}

#' Calculate index summaries
#' 
#' Generates some simple summaries for each unique index.
#' 
#' @import tidyverse
#' @param df A dataframe which must contain the column "n" which represents the count of a read.
#' @param index_col The name of the column contain the index barcodes as a string. This column must be present in "df".
#' @param valid_indices. A vector of all the valid indices for "index_col".
#' @returns A dataframe with the follow columns:
#'         - index_col: String, The column containing the index barcodes.
#'         - idx_n: Numeric, Number of reads associated with a specific index barcode.
#'         - fraction: Numeric, "idx_n" divided by the total number of reads in the run.
#'         - expected: Boolean, True if the index barcode is in "valid_indices" otherwise False.
#'         - contains_n: Boolean, True if the index barcode contains "N" in its sequence, otherwise False.
#'         - lv_dist: Numeric, Edit distance from a valid index barcode.
#'         - ham_dist: Numeric, Hamming distance from a valid index barcode.
get_index_summary= function(df, index_col, valid_indices) {
  output_summary= df %>% dplyr::group_by(pick(all_of(index_col))) %>% 
    dplyr::summarise(idx_n= sum(n)) %>% dplyr::ungroup() %>%
    dplyr::mutate(fraction= round(idx_n/sum(idx_n), 5),
                  expected= ifelse(.[[index_col]] %chin% valid_indices, T, F),
                  contains_n= ifelse(grepl('N', .[[index_col]]), T, F),
                  lv_dist= apply(stringdist::stringdistmatrix(.[[index_col]], valid_indices, method="lv"), 1, min),
                  ham_dist= apply(stringdist::stringdistmatrix(.[[index_col]], valid_indices, method="hamming"), 1, min)) %>%
    dplyr::arrange(desc(fraction))
  return(output_summary)
}

#' Total counts barplot
#' 
#' Creates the total counts barplot with bars colored by the barcode type,
#' either a cell line barcode or control barcode.
#'
#' @import tidyverse
#' @param filtered_counts Filtered counts dataframe.
#' @param id_cols Vector of columns names that identify each sample.
#' @param facet_col String name of the column in filtered_counts to facet the plot.
#'                  This can be left as NA if there isn't a column to facet on.
#' @returns Returns a ggplot object.
create_total_counts_barplot= function(filtered_counts, id_cols, facet_col= NA) {
  # Validation: Check that id_cols and facet_col exist in filtered counts.
  if(!validate_columns_exist(na.omit(c(id_cols, facet_col)), filtered_counts)) {
    stop('Some input columns were not detected in filtered counts.')
  }
  
  # Sum up reads 
  total_counts= filtered_counts %>%
    dplyr::mutate(barcode_type = case_when(
    !("cb_name" %in% names(cur_data())) ~ "cell line", is.na(cb_name) ~ "cell line", TRUE ~ "ctrl barcode")) %>%
    tidyr::unite(all_of(id_cols), col= 'sample_id', sep= ':', remove= FALSE, na.rm= FALSE) %>%
    dplyr::group_by(pick(all_of(na.omit(c('sample_id', facet_col, 'barcode_type'))))) %>%
    dplyr::summarise(total_counts= sum(n)) %>% dplyr::ungroup()
  
  # Create total counts plot
  total_counts_plot= total_counts %>% 
    ggplot(aes(x= sample_id, y= total_counts, fill= barcode_type)) +
    geom_col(alpha= 0.75, position= 'identity') +
    geom_hline(yintercept= 10^4, linetype= 2) + 
    {if(!is.na(facet_col)) facet_wrap(~.data[[facet_col]], scale= 'free_x')} +
    labs(x= 'Sample constructed using id_cols', y= 'Total counts', fill= 'Barcode\ntype', 
         title= 'Filtered counts - unstacked') + 
    theme_bw() + theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(total_counts_plot)
}

#' Cell line recover barplot
#' 
#' Creates barplots of the cell lines recovered. The parameter "plot_type" can be used to plot the percentage or
#' the total cell line counts on teh y axis. The parameter "include_ctrl_bcs" can be used to include the control
#' barcodes in the cell line count. 
#' 
#' @import tidyverse
#' @param filtered_counts Filtered counts dataframe.
#' @param id_cols Vector of column names that identify each sample.
#' @param facet_col String name of the column in filtered_counts to facet the plot.
#' @param value_col String name of the column in filtered_counts that contains the counts.
#' @param count_threshold Threshold used to determine low counts.
#' @param plot_type String of either "percent" or "count" to adjust the y axis to be either the percentage or the 
#'                  total number of cell lines.
#' @param include_ctrl_bcs Boolean. Set to TRUE if control barcodes are to be counted. 
#' @returns Returns a ggplot plot.
create_recovery_barplot= function(filtered_counts, id_cols, facet_col= NA, value_col= 'n', count_threshold, 
                                  plot_type= 'percent', include_ctrl_bcs= FALSE) {
  # Validation: Check that id_cols, facet_col, or value_col exist in filtered counts.
  if(!validate_columns_exist(na.omit(c(id_cols, facet_col, value_col)), filtered_counts)) {
    stop('In create_recovery_barplot, some required input columns were not detected.')
  }
  
  # Filter out control barcodes if it is specified.
  if(include_ctrl_bcs == FALSE) {
    filtered_counts= filtered_counts %>% filter_control_barcodes()
  }
  
  # Count number of cell lines/ barcodes for a detection group.
  recovery= filtered_counts %>%
    dplyr::add_count(pick(all_of(id_cols)), name= 'total_num_cls') %>%
    dplyr::mutate(detect_type= case_when(.data[[value_col]] == 0 ~ 'Not detected',
                                         .data[[value_col]] <= count_threshold ~ 'Low counts',
                                         .data[[value_col]] > count_threshold ~ 'Detected')) %>%
    dplyr::count(pick(all_of(na.omit(c(id_cols, facet_col, 'detect_type', 'total_num_cls')))), name= 'num_cls_by_type') %>%
    tidyr::unite(all_of(id_cols), col= 'sample_id', sep= ':', remove= FALSE, na.rm= FALSE) %>%
    dplyr::mutate(percent= (num_cls_by_type / total_num_cls) * 100)
  
  # Set the y axis depending on the plot type.
  if(plot_type == 'count') {
    y_col= 'num_cls_by_type'
    y_text= 'Number of cell lines'
  } else {
    if(plot_type != 'percent') {
      print(paste0('Warning: ', plot_type, ' is not a valid plot type. Please use either count or percent.'))
      print('Defaulting to percent plot.')
    }
    y_col= 'percent'
    y_text= 'Percentage of cell lines recovered (%)'
  }
  
  # Create recovery plot.
  recov_plot= recovery %>%
    ggplot(aes(x= sample_id, y= .data[[y_col]], fill= reorder(detect_type, dplyr::desc(detect_type)))) +
    geom_col(position= 'stack', alpha= 0.75) +
    {if(!is.na(facet_col)) facet_wrap(~.data[[facet_col]], scale= 'free_x')} +
    labs(x= 'Sample constructed using id_cols', y= y_text, fill= '', title= 'Cell line recovery') + 
    theme_bw() + theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(recov_plot)
}

#' Cumulative reads plot
#' 
#' Creates a line plot of the cumulative reads.
#' 
#' @import tidyverse
#' @param input_df Input dataframe. Usually is the filtered_counts dataframe.
#' @param id_cols Vector of column names that identify every PCR well.
#' @param counts_col Name of the column that contains the values. Defaults to "n".
#' @param mark1 Percentage of reads to mark. Draws a line at a specified percentage to indicate the number of 
#'              cell lines needed to reach this percentage of reads. Defaults to 0.5.
#' @param mark2 Percentage of reads to mark. Draws a line at a specified percentage to indicate the number of 
#'              cell lines needed to reach this percentage of reads. This parameter should be greater than the 
#'              value specified for "mark1". Defaults to 0.95.
#' @param contains_cb Boolean. If control barcodes are used, this can be set to TRUE so that points 
#'                    corresponding to the control barcodes will be colored on the plot. Defaults to FALSE.
#' @param order_aucs Boolean, when there are multiple facets, this can be set to TRUE to sort the facets by 
#'                   the auc value. The aucs will be sorted in descending order. Defaults to FALSE.
#' @returns Returns a ggplot object.
create_cdf_plot= function(input_df, id_cols, counts_col= 'n', mark1= 0.5, mark2= 0.95, 
                          contains_cbs= FALSE, order_aucs= FALSE) {
  # Validation: Check that id_cols and counts_col are in the input dataframe.
  if(!validate_columns_exist(c(id_cols, counts_col), input_df)) {
    stop('In create_cdf_plot, some required input columns were not detected.')
  }
  
  # Validation: mark1 should be less than mark2.
  if(mark1 > mark2 | mark1 < 0 | mark1 > 1) {
    stop('Mark values must be between 0 and 1. mark1 should be less than mark2')
  }
  
  # Determine percentages, ranks and cumulative percentages
  calc_cummulative= input_df %>% dplyr::group_by(pick(all_of(id_cols))) %>%
    dplyr::arrange(dplyr::desc(.data[[counts_col]])) %>%
    dplyr::mutate(expected_num_cls= dplyr::n(),
                  total_counts= sum(.data[[counts_col]]), pct_counts= .data[[counts_col]]/total_counts,
                  cum_pct= cumsum(pct_counts), 
                  rank= row_number(), rank_pct= rank / expected_num_cls) %>% dplyr::ungroup()
  
  # Find the number of cell lines needed to reach mark1 and mark2
  mark1_values= calc_cummulative %>% dplyr::filter(cum_pct >= mark1) %>% 
    dplyr::group_by(pick(all_of(id_cols))) %>% dplyr::arrange(cum_pct) %>% 
    dplyr::filter(row_number() == 1) %>% dplyr::ungroup() %>% 
    dplyr::select(all_of(id_cols), rank_pct= rank_pct, mark1_rank= rank, mark1_loc= rank_pct)
  mark2_values= calc_cummulative %>% dplyr::group_by(pick(all_of(id_cols))) %>% 
    dplyr::mutate(auc= sum(cum_pct * (1 / expected_num_cls))) %>% # calculate aucs
    dplyr::filter(cum_pct >= mark2) %>% dplyr::arrange(cum_pct) %>%
    dplyr::filter(row_number() == 1) %>% dplyr::ungroup() %>% 
    dplyr::select(all_of(id_cols), rank_pct= rank_pct, mark2_rank= rank, mark2_loc= rank_pct, auc)
  
  # Create cdf plot
  data_for_plot= calc_cummulative %>% 
    dplyr::left_join(mark1_values, by= c(id_cols, 'rank_pct')) %>% 
    dplyr::left_join(mark2_values, by= c(id_cols, 'rank_pct')) %>%
    tidyr::unite(all_of(id_cols), col= 'facet_name', sep= ':', remove= TRUE, na.rm= FALSE)
  
  # Reorder by aucs if specified
  if(order_aucs) {
    data_for_plot= data_for_plot %>% dplyr::arrange(dplyr::desc(auc)) %>%
      dplyr::mutate(facet_name= base::factor(facet_name, levels= unique(facet_name)))
  }
  
  # Create plot
  output_plot= data_for_plot %>%
    ggplot(aes(x= rank_pct, y= cum_pct)) +
    # Color control barcodes if specified
    {if(contains_cbs) geom_point(. %>% dplyr::filter(!is.na(cb_name)),
                                 mapping= aes(x= rank_pct, y=cum_pct, color= reorder(cb_name, cb_log2_dose)), size= 2)} + 
    geom_line(color='black') +
    # point for mark1 of counts
    geom_segment(aes(x= -Inf , y= mark1, xend= mark1_loc, yend = mark1), color= 'black', linetype= 2) +
    geom_segment(aes(x= mark1_loc, y= -Inf, xend = mark1_loc, yend = mark1), color= 'black', linetype= 2) +
    geom_label(aes(x= mark1_loc, y= 0.25, label= mark1_rank), hjust= 0, color= 'black') +
    # point for 95% of counts
    geom_segment(aes(x= -Inf , y= mark2, xend= mark2_loc, yend= mark2), color= 'black', linetype= 2) +
    geom_segment(aes(x= mark2_loc, y= -Inf, xend= mark2_loc, yend= mark2), color= 'black', linetype= 2) +
    geom_label(aes(x= mark2_loc, y= 0.75, label= mark2_rank), hjust= 0, color= 'black') +
    # label for auc
    #geom_label(aes(x= mark2_loc, y= 0.1, label= paste0('auc ', round(auc, 3))), hjust= 'inward', color= 'black') +
    geom_label(. %>% dplyr::filter(!is.na(auc)), mapping= aes(label= paste0('auc ', round(auc, 3))),
               x= 1, y= 0.25, hjust= 'inward', vjust= 'inward', color= 'black') +
    facet_wrap(~facet_name) + 
    labs(x= '% rank of unique reads', y= 'Cumulative percentage', color= 'CBs') + theme_bw()
  
  return(output_plot)
}

#' Control barcode scatter plot
#' 
#' Creates a scatter plot of the control barcodes.
#' 
#' @import tidyverse
#' @param normalized_counts Dataframe output from the normalize module.
#' @param id_cols Vector of column names that identify every PCR well.
#' @param value_col Name of the column that contains the values.
#' @param pseudocount a pseudo raw count used for missing cell line counts, also used when computing log2 on raw counts
#' @returns Returns a ggplot object.
create_ctrlBC_scatterplots= function(normalized_counts, id_cols, value_col= 'log2_n', pseudocount = 20) {
  # Validation: Check that id_cols and value_col exist in filtered counts.
  if(value_col == "log2_n" & validate_columns_exist(c(id_cols, "n"), normalized_counts) & 
     !(validate_columns_exist(c("log2_n"), normalized_counts))){
    print("Calculating log2_n")
    normalized_counts %<>% mutate(log2_n = log2(n+1))
  }else if(!validate_columns_exist(c(id_cols, value_col), normalized_counts)) {
    stop('Some input columns were not detected in normalized counts.')
  }
  
  # Detect norm_r2 and norm_mae. If columns do not exist, then roughly calculate those columns.
  if(any(!c('norm_r2', 'norm_mae') %in% colnames(normalized_counts))) {
    print('WARNING: Columns "norm_r2" and/or "norm_mae" were not detected in normalized_counts.', quote= FALSE)
    print('Calculating both columns - this method may not be as robust as the normalize module.')
    
    normalized_counts= normalized_counts %>%
      filter_control_barcodes() %>%
      dplyr::filter(!is.na(cb_ladder), n != 0) %>%
      dplyr::group_by(pick(all_of(id_cols))) %>%
      dplyr::mutate(mean_y= mean(cb_log2_dose),
                    residual2= (cb_log2_dose - log2_n)^2,
                    squares2= (cb_log2_dose - mean_y)^2,
                    norm_r2= 1 - sum(residual2) / sum(squares2),
                    norm_mae= median(abs(cb_log2_dose - log2_n))) %>% ungroup()
  } 
  
  # Filter for just the control barcodes, create a profile_id for faceting, 
  # and determine the x and y positions for the r2 + mae label.
  cb_trend= normalized_counts %>% dplyr::filter(!is.na(cb_name), !is.na(cb_ladder)) %>%
    tidyr::unite(all_of(id_cols), col= 'profile_id', sep= ':', remove= TRUE) %>%
    dplyr::group_by(profile_id) %>% dplyr::mutate(label_x_pos= min(.data[[value_col]]),
                                                  label_y_pos= max(cb_log2_dose)) %>% dplyr::ungroup()
  
  # Create control barcode trend plot
  trend_scatter_plot= cb_trend %>% ggplot(aes(x= .data[[value_col]], y= cb_log2_dose)) + 
    geom_point() +
    geom_abline(aes(slope=1, intercept= cb_intercept) , color='blue', alpha= 0.5) +
    geom_text(aes(x= label_x_pos, y= label_y_pos,
                  label= paste0('r2= ', round(norm_r2, 4), '\nmae= ', round(norm_mae, 4))), 
              hjust='inward', vjust='inward', alpha= 0.5) +
    facet_wrap(~profile_id, scales='free_x') +
    labs(title= 'Linear fit of control barcodes across all samples') + theme_bw()
  
  return(trend_scatter_plot)
}

#' Heatmap of correlations
#' 
#' Creates a correlation heatmap. A matrix of values is created from the input_df. The row_id_cols
#' are used identify each row and the col_id_cols are used to identify each column. The value_col is 
#' used to fill the matrix. Correlations are then computed.
#' 
#' @import tidyverse
#' @import WGCNA
#' @import reshape2
#' @import scales
#' @param input_df Dataframe.
#' @param row_id_cols Vector of column names from input_df that identifies the cell lines. For example,
#'                    this can be "depmap_id", "ccle_name" if only cell lines exist. It can also be 
#'                    "depmap_id", "ccle_name", "cb_name" if control barcodes are also present.
#' @param col_id_cols Vector of column names from input_df that identifies the PCR wells or conditions.
#'                    For example, this can be "pcr_plate", "pcr_well" or a list of conditions like those in sig_cols.
#' @param value_col String name of the column in input_df to be used as the values.
#' @param cor_method WGCNA correlation method. This defaults to "pearson".
#' @returns Returns a ggplot object.
create_cor_heatmap= function(input_df, row_id_cols, col_id_cols, value_col,
                             cor_method= 'pearson') {
  
  # Validate that specified columns are in the dataframe.
  if(!validate_columns_exist(c(row_id_cols, col_id_cols, value_col), input_df)) {
    stop('Not all columns were detected in the input dataframe.')
  }
  
  # Create row and column names for pivoting to a matrix
  correlation_mx= input_df %>%
    tidyr::unite(all_of(row_id_cols), col= 'row_id', sep= ':', remove= TRUE) %>%
    tidyr::unite(all_of(col_id_cols), col= 'col_id', sep= ':', remove= TRUE)
  
  # Check that the row and column ids specify one value.
  validate_ids= correlation_mx %>% dplyr::group_by(row_id, col_id) %>%
    dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
  if(nrow(validate_ids) != 0) {
    print('The provide columns specify more than one value.')
    print(head(validate_ids))
    stop('Multiple values detected for a unique combination of "row_id_cols" and "col_id_cols".')
  }
  
  # Pivot and calculate correlations
  correlation_mx= correlation_mx %>% reshape2::acast(row_id~col_id, value.var= value_col) %>%
    WGCNA::cor(use= 'pairwise.complete.obs', method= cor_method)
  
  # Create heatmap 
  cor_heatmap= correlation_mx %>% reshape2::melt() %>% 
    ggplot(aes(x= Var1, y= Var2, fill= value)) + 
    geom_tile() +
    scale_fill_gradientn(breaks= c(0, 0.5, 1), 
                         colours= c('blue', 'white','red'),
                         limits=c(0, 1), oob= scales::squish) +
    labs(x= '', y= '', fill= '', title= paste0('Correlations using ', value_col)) +
    theme_bw() + theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(cor_heatmap)
}

#' Scatter plots of two replicates
#' 
#' From a long table, creates scatter plots to two replicates.
#' 
#' @import tidyverse
#' @import ggpmisc
#' @param input_df Dataframe.
#' @param cell_line_cols List of column names used to identify each cell line or control barcode.
#' @param replicate_group_cols List of column names that describe a group of similar conditions.
#' @param replicate_col Name of the column that specifies the replicate. This column should not be 
#'                      in replicate_group_cols!
#' @param value_col Name of the column in input_df that contains the values.
#' @param x_axis_rep String of the replicate identifier that should be on the x axis of the plot.
#' @param y_axis_rep String of the replicate identifier that should be on the y axis of the plot.
#' @returns Returns a ggplot object or NULL if all entries are filtered out.
create_replicate_scatterplots= function(input_df, cell_line_cols, replicate_group_cols, replicate_col, value_col,
                                        x_axis_rep= '1', y_axis_rep= '2') {
  # Validation: Check that input columns are present in the dataframe.
  if(!validate_columns_exist(c(cell_line_cols, replicate_group_cols, replicate_col, value_col), input_df)) {
    stop('Some input columns were not detected in normalized counts.')
  }
  
  # Validation: Check that replicate_col is not in replicate_group_cols.
  if(replicate_col %in% replicate_group_cols) {
    stop(paste0(replicate_col, ' should not be included in replicate_group_cols!'))
  }
  
  reps_piv= input_df %>% 
    tidyr::unite(all_of(replicate_group_cols), col= 'replicate_group', sep= ':', remove= TRUE, na.rm= FALSE) %>%
    dplyr::group_by(pick(all_of(c(cell_line_cols, 'replicate_group')))) %>%
    dplyr::filter(!is.na(.data[[replicate_col]]), .data[[replicate_col]] != '', .data[[value_col]] != 0,
                  dplyr::n() >= 2) %>% 
    dplyr::ungroup() 
  
  # Retun a null object if no entries pass the filter.
  if(nrow(reps_piv) == 0) {return(NULL)}
  
  reps_piv= reps_piv %>%
    pivot_wider(id_cols= all_of(c(cell_line_cols, 'replicate_group')),
                names_from= replicate_col, names_prefix= replicate_col, values_from= value_col) %>%
    dplyr::mutate(type= ifelse(!is.na(lua), "cell line", "control barcode")) %>% dplyr::ungroup()
  
  # Create names of the columns to plot on xy axes
  x_col_name= paste0(replicate_col, x_axis_rep)
  y_col_name= paste0(replicate_col, y_axis_rep)
  
  reps_scatter= reps_piv %>% dplyr::filter(!is.na(.data[[x_col_name]]), !is.na(.data[[y_col_name]])) %>%
    ggplot(aes(x= .data[[x_col_name]], y= .data[[y_col_name]])) +
    geom_abline(color='black', linewidth=0.5, linetype=2) +
    geom_point(aes(color= type), alpha=0.75) +
    ggpmisc::stat_correlation(mapping = use_label(c("R2", "n")))+ 
    facet_wrap(~replicate_group, scales= 'free') +
    labs(x= paste0(replicate_col, '1 ', value_col), y= paste0(replicate_col, '2 ', value_col),
         title= paste0('Scatter plot of ', replicate_col, ' with ', value_col)) +
    theme_bw()
  
  return(reps_scatter)
}

#' QC_images
#'
#' Takes in various pipeline outputs and generates 11 QC files.
#'
#' @param raw_counts_uncollapsed_path Path to the raw_counts_uncollapsed file.
#' @param prism_barcode_counts Dataframe of prism barcodes identified in the run.
#' @param unknown_barcode_counts Dataframe of unknown barcodes.
#' @param annotated_counts Annotated counts dataframe from the filter_raw_reads module.
#' @param normalized_counts Normalized counts dataframe from the normalize module. This is an optional parameter.
#' @param l2fc L2FC dataframe from the compute_l2fc module. This is used for the bio_reps plot. 
#' @param sample_meta Dataframe of the sample metadata for the sequencing run.
#' @param barcode_col String name of the column containing the barcode sequences.
#' @param cell_line_cols Vector of sample meta column names used to describe a cell line or barcode.
#' @param id_cols Vector of sample meta column names used to identify each PCR well. 
#'                This defaults to "pcr_plate", "pcr_well".
#' @param sig_cols Vector of sample meta column names used to identify a unique treatment condition.
#' @param control_type String of how the negative controls are designated in the pert_type column in the sample_meta.
#' @param count_threshold Threshold for low read counts.
#' @param reverse_index2 Boolean set to TRUE if the sequencing involved the reverse complement workflow.
#' @param out Path to the directory to save the QC images.
#' @param pseudocount a pseudo raw count used for missing cell line counts, also used when computing log2 on raw counts
#' @returns NA. QC images are written out to the specified folder.
QC_images= function(raw_counts_uncollapsed_path,
                    prism_barcode_counts, unknown_barcode_counts,
                    annotated_counts, normalized_counts= NA, l2fc, 
                    sample_meta,
                    barcode_col= 'forward_read_barcode',
                    id_cols= c('pcr_plate', 'pcr_well'),
                    cell_line_cols= c('depmap_id', 'lua'), 
                    sig_cols,
                    control_type= 'ctl_vehicle', count_threshold= 40, 
                    chunk_size= 10^6, pseudocount = 20,
                    reverse_index2= FALSE, out= NA) {
  
  # Required packages ----
  require(tidyverse)
  require(magrittr)
  require(data.table)
  require(reshape2)
  require(WGCNA)
  require(scales)
  require(ggpmisc)
  
  # Some preprocessing ----
  # Set out directory if none is specified.
  if(is.na(out)) {out= getwd()}
  
  # Create empty vector to collect potential errors when running QCs
  skipped_qcs= c()
  
  # Count number of distinct profile to help scale some plots.
  num_profiles= annotated_counts %>% dplyr::distinct(pick(all_of(id_cols))) %>% nrow()
  
  # Detect if there are wells with control barcodes.
  cb_check= sample_meta %>% dplyr::filter(!cb_ladder %in% c('None', 'empty', '', ' ', NA),
                                          !pert_type %in% c("empty", "", "CB_only", NA))
  contains_cbs= nrow(cb_check) > 0 # set boolean
  
  # Create filtered_counts df from annotated_counts 
  filtered_counts= annotated_counts %>% dplyr::filter(expected_read)
  
  # Sequencing QCs ____________________ ----
  ## 1. Purity metrics ----
  print('1. Generating QC table ...')
  create_qc_table(raw_counts_uncollapsed_path= raw_counts_uncollapsed_path, 
                  unknown_barcode_counts= unknown_barcode_counts,
                  prism_barcode_counts= prism_barcode_counts,
                  filtered_counts= filtered_counts,
                  value_col= 'n', chunk_size= chunk_size,
                  output_path= paste0(out, '/QC_table.csv'))
  
  ## 2. Index count summaries ----
  print('2. Generating index counts tables ...')

  # Pull out headers to perform checks
  raw_counts_uncollapsed_headers= read_data_table(raw_counts_uncollapsed_path, nrows= 0)
  
  # Check that "index_1" is present. If so, calculate index summary and write out.
  if('index_1' %in% colnames(sample_meta) & 'index_1' %in% colnames(raw_counts_uncollapsed_headers)) {
    # Aggregate over index_1 using chunks
    # Action is set to a data.table summarize with summing
    index1_chunks= process_in_chunks(large_file_path= raw_counts_uncollapsed_path, chunk_size= chunk_size, 
                                     action= function(x) x[, list(n= sum(n)), by= index_1])
    
    # Create vector of unique index_1 values
    expected_index1= unique(sample_meta$index_1)
    
    # Call get_index_summary over index1_chunks as a full table, then write out table
    index1_counts= get_index_summary(data.table::rbindlist(index1_chunks), 'index_1', expected_index1)
    index1_counts %>% write.csv(file= paste(out, 'index1_counts.csv', sep= '/'), row.names= FALSE, quote= FALSE)
  } else {
    print('Column "index_1" not detected. Skipping index 1 summaries ...', quote= FALSE)
  }
  
  # Do the same for index 2.
  # Reverse index 2 barcodes if it is indicated and if "index_2" exists
  if(reverse_index2 & 'index_2' %in% colnames(sample_meta) ) {
    print("Reverse-complementing index 2 barcode.")
    sample_meta$index_2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
  }
  
  if('index_2' %in% colnames(sample_meta) & 'index_2' %in% colnames(raw_counts_uncollapsed_headers)) {
    # Aggregate over index_2 using chunks
    # Action is set to a data.table summarize with summing
    index2_chunks= process_in_chunks(large_file_path= raw_counts_uncollapsed_path, chunk_size= chunk_size, 
                                     action= function(x) x[, list(n= sum(n)), by= index_2]) 
    
    # Create vector of unique index_2 values
    expected_index2= unique(sample_meta$index_2)
    
    # Call get_index_summary over index2_chunks as a full table, then write out table
    index2_counts= get_index_summary(data.table::rbindlist(index2_chunks), 'index_2', expected_index2)
    index2_counts %>% write.csv(file= paste(out, 'index2_counts.csv', sep= '/'), row.names= FALSE, quote= FALSE)
  } else {
    print('Column "index_2" not detected. Skipping index 2 summaries ...',  quote= FALSE)
  }
  
  ## 3. Total counts ----
  print('3. Generating total_counts image ...')
  potential_error= base::tryCatch({
    tc= create_total_counts_barplot(filtered_counts, id_cols, facet_col= 'pcr_plate')
    
    pdf(file=paste(out, "total_counts.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
    print(tc)
    dev.off()
    rm(tc)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the total counts barplot. Skipping this output ...') 
    return('Totalc ounts image')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## 4. Unknown barcodes ----
  print('4. Generating table of unknown barcode reads ...')
  potential_error= base::tryCatch({
    unknown_totals= unknown_barcode_counts[, .(well_total= sum(n)), by= id_cols]
    prism_totals= prism_barcode_counts[, .(well_total= sum(n)), by= id_cols]
    well_totals= data.table::rbindlist(list(unknown_totals, prism_totals))[, .(well_total= sum(well_total)), by= id_cols]
    
    unknown_barcodes= unknown_barcode_counts %>% 
      dplyr::filter(.data[[barcode_col]] != 'unknown_reads') %>%
      dplyr::left_join(well_totals, by= id_cols) %>%
      dplyr::mutate(read_percent= n / well_total) %>%
      dplyr::group_by(pick(all_of(barcode_col))) %>%
      dplyr::summarise(total= sum(n),
                       median_read= median(n),
                       median_percent= median(read_percent),
                       max_read= max(n),
                       max_percent= max(read_percent),
                       num_wells= dplyr::n()) %>% dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(median_percent))
    
    unknown_barcodes %>% write.csv(file= paste(out, 'unknown_barcodes_summary.csv', sep= '/'), 
                                   row.names= FALSE, quote= FALSE)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the summary unknown barcode reads. Skipping this output ...') 
    return('unknown barcode reads')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  #
  
  # Assay QCs _________________________ ----
  ## 5. Cell lines recovered ----
  print('5. Generating cell_lines_present image ...')
  potential_error= base::tryCatch({
    cl_rec= create_recovery_barplot(filtered_counts, id_cols= id_cols, facet_col= 'pcr_plate', 
                                    count_threshold= count_threshold, plot_type= 'percent')
    
    pdf(file= paste(out, "cell_lines_present.pdf", sep="/"),
        width= sqrt(num_profiles)*2, height= sqrt(num_profiles))
    print(cl_rec)
    dev.off()
    rm(cl_rec)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the recovery barplot. Skipping this output ...')
    return('CL recovery')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## 6. Cell line contaminants ----
  print('6. Generating cell line contaminants ...')
  potential_error= base::tryCatch({
    contams= annotated_counts %>% dplyr::filter(expected_read == FALSE) %>%
      dplyr::mutate(barcode_id= ifelse(is.na(depmap_id), cb_name, depmap_id)) %>%
      dplyr::group_by(forward_read_barcode, barcode_id) %>%
      dplyr::summarise(num_wells= n(), median_n= median(n), max_n= max(n)) %>% ungroup() %>%
      dplyr::arrange(desc(num_wells))
    
    contams %>% write.csv(file= paste(out, 'contam_cell_lines.csv', sep='/'), row.names= FALSE, quote= FALSE)
    rm(contams)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the contaminants file. Skipping this output ...') 
    return('contam_cell_lines.csv')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## 7. Cumulative counts by lines in negcons ----
  print('7. Generating cumulative image ...')
  potential_error= base::tryCatch({
    cdf_plot= create_cdf_plot(filtered_counts %>% dplyr::filter(pert_type == control_type), 
                              id_cols= id_cols, 
                              counts_col= 'n', 
                              mark1= 0.5, mark2= 0.95, 
                              contains_cbs= contains_cbs, order_aucs= TRUE) +
      labs(title= 'Cumulative reads in negative controls.')
    
    pdf(file=paste(out, 'cdf_plot.pdf', sep= '/'),
        width= sqrt(num_profiles) * 2, height= sqrt(num_profiles))
    print(cdf_plot)
    dev.off()
    rm(cdf_plot)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the cdf plot. Skipping this output ...') 
    return('cdf plot')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs= c(skipped_qcs, potential_error)
  }
  
  ## 8. Control barcode trends ----
  if(contains_cbs & is.data.frame(normalized_counts)) {
    print('8. Generating control_barcode_trend image')
    potential_error= base::tryCatch({
      trend_sc= create_ctrlBC_scatterplots(normalized_counts %>% dplyr::filter(!cb_ladder %in% c(NA, FALSE, 'none', '')), 
                                           id_cols, value_col= 'log2_n',
                                           pseudocount = pseudocount)
      
      pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
          width= sqrt(num_profiles) * 2, height= sqrt(num_profiles) * 2)
      print(trend_sc)
      dev.off()
      rm(cb_trend, trend_sc)
    }, error= function(e) {
      print(e)
      print('Encountered an error when creating the CB trends plot. Skipping this output ...') 
      return('cb trend')
    })
    
    # Collect returned string if an error occurred
    if(!is.null(potential_error)) {
      skipped_qcs= c(skipped_qcs, potential_error)
    }
  } else {
    print('8. No control barcodes detected. Skipping control_barcode_trend image.')
  }
  
  ## 9. Sample correlation -----
  print('9. Generating sample_cor image ...')
  potential_error= base::tryCatch({
    cor_df= filtered_counts %>%
      filter_control_barcodes() %>%
      dplyr::filter(!pert_type %in% c(NA, 'empty', '', 'CB_only')) %>%
      dplyr::mutate(log2_n= log2(n + 1))
    cp= create_cor_heatmap(input_df= cor_df,
                           row_id_cols= c('depmap_id'),
                           col_id_cols= c(sig_cols, id_cols),
                           value_col= 'log2_n')
    
    pdf(file= paste(out, 'sample_cor.pdf', sep= '/'),
        width= sqrt(num_profiles) * 2, height= sqrt(num_profiles) * 2)
    print(cp)
    dev.off()
    rm(correlation_matrix, cp)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the sample_cor figure. Skipping this output ...')
    return('sample_cor')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs= c(skipped_qcs, potential_error)
  }
  
  ## 10. Tech rep correlations ----
  if(is.data.frame(normalized_counts) & 'tech_rep' %in% colnames(normalized_counts)) {
    # Check if there are more at least two tech reps
    unique_tech_reps= na.omit(unique(normalized_counts$tech_rep))
    
    if(length(unique_tech_reps) >= 2) {
      print('10. Generating tech rep correlations image ...')
      # Set up replicate groups depending "bio_rep" column
      if('bio_rep' %in% colnames(normalized_counts) & !'bio_rep' %in% sig_cols) {
        replicate_group_cols= c(sig_cols, 'bio_rep')
      } else {
        replicate_group_cols= sig_cols
      }
      
      # Handle cases if control barcodes are used.
      if('cb_name' %in% colnames(normalized_counts)) {
        unique_cell_line_cols= c(cell_line_cols, 'cb_name')
      } else {
        unique_cell_line_cols= cell_line_cols
      }
      
      # Create replicate scatter plot
      potential_error= base::tryCatch({
        tech_reps_plt= create_replicate_scatterplots(input_df= normalized_counts, 
                                                     cell_line_cols= unique_cell_line_cols, 
                                                     replicate_group_cols= replicate_group_cols, 
                                                     replicate_col= 'tech_rep', 
                                                     value_col= 'log2_n')
        if(!is.null(tech_reps_plt)) {
          pdf(file=paste(out, "tech_reps_plt.pdf", sep="/"),
              width=sqrt(num_profiles), height=sqrt(num_profiles))
          print(tech_reps_plt)
          dev.off()
        } else {
          print('No technical replicates detected - skipping plot.')
        }
      }, error= function(e) {
        print(e)
        print('Encountered an error when creating the tech_reps_plt figure. Skipping this output ...')
        return('tech_reps_plt')
      })
      
      # Collect returned string if an error occurred
      if(!is.null(potential_error)) {
        skipped_qcs= c(skipped_qcs, potential_error)
      }
      
    } else {
      print('10. No technical replicates detected. Skipping tech_reps scatter plot.')
    }
  } else {
    print('10. No technical replicates detected. Skipping tech_reps scatter plot.')
  }
  
  ## 11. Bio rep correlations ----
  if('bio_rep' %in% colnames(l2fc)) {
    unique_bio_reps= na.omit(unique(l2fc$bio_rep))
    
    if(length(unique_bio_reps) >= 2) {
      l2fc_with_log2= l2fc %>% dplyr::mutate(log2_mean_normalized_n= log2(mean_normalized_n))
      
      # Bio replicate scatter plots
      # This is just another visualization that isn't being used.
      # bio_reps_plt= create_replicate_scatterplots(input_df= l2fc_with_log2s, 
      #                                           cell_line_cols= cell_line_cols, 
      #                                           replicate_group_cols= sig_cols, 
      #                                           replicate_col= 'bio_rep', 
      #                                           value_col= 'log2_mean_normalized_n')
      # if(!is.null(bio_reps_plt)) {
      #   pdf(file=paste(out, "bio_reps_plt.pdf", sep="/"),
      #       width=sqrt(num_profiles), height=sqrt(num_profiles))
      #   print(bio_reps_plt)
      #   dev.off()
      # } else {
      #   print('No technical replicates detected - skipping plot.')
      # }
      
      # Bio replicate heatmap
      print('11. Generating bio rep correlations heatmap ...')
      potential_error= base::tryCatch({
        bio_corr_hm= create_cor_heatmap(input_df= l2fc_with_log2, 
                                        row_id_cols= cell_line_cols, 
                                        col_id_cols= c(sig_cols, 'bio_rep'), 
                                        value_col= 'l2fc',
                                        cor_method= 'pearson') 
        pdf(file= paste(out, 'bio_corr_hm.pdf', sep= '/'),
            width= sqrt(num_profiles), height= sqrt(num_profiles))
        print(bio_corr_hm)
        dev.off()
      }, error= function(e) {
        print(e)
        print('Encountered an error when creating the bio_corr_hm figure. Skipping this output ...') 
        return('bio_corr_hm')
      })
      
      # Collect returned string if an error occurred
      if(!is.null(potential_error)) {
        skipped_qcs= c(skipped_qcs, potential_error)
      }
      
    } else {
      print('11. No biological replicates detected. Skipping bio_rep heatmap.')
    }
  }
  
  # End _________________________ ----
  print('QCs finishing!')
  if(length(na.omit(skipped_qcs)) != 0) {
    print(paste0('WARNING: The following ', length(skipped_qcs), ' QCs encountered errors and were skipped - '))
    print(na.omit(skipped_qcs))
  } else {
    print('No errors encountered.')
  }
}
