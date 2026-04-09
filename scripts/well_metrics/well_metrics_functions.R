options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(PRROC)
library(dplyr)

#' Compute read stats
#'
#' This function calculates summary statistics related to reads and cell line recovery for annotated counts data,
#' combining metadata and performing group-wise computations.
#'
#' @param unknown_counts A data.table containing reads that did not map to any PRISM barcodes.
#' @param annotated_counts A data frame containing read data with annotations for barcodes, cell lines, and counts.
#' @param id_cols A character vector of column names that define a unique PCR well.
#'   Default: `c("pcr_plate", "pcr_well")`.
#' @param cell_line_cols A character vector of column names that define a unique cell line.
#'   Default: `c("cell_set", "pool_id", "depmap_id", "lua")`.
#' @param metadata_cols A character vector of column names of additional coluums to include.
#'   Default: 'c("pert_type", "pert_plate", "cell_set")'.
#' @param metric A string of the column name of the metric to use for calculations. Default: `"n"`.
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered
#'   Default: `40`.
#'
#' @return A data frame summarizing read statistics for each group, including total reads, expected reads,
#'         control barcode reads, recovered cell lines, and their fractions.
compute_read_stats = function(unknown_counts, annotated_counts,
                              id_cols = c("pcr_plate", "pcr_plate"),
                              cell_line_cols = c("cell_set", "pool_id", "depmap_id", "lua"),
                              metadata_cols = c("pert_type", "pert_plate", "cell_set"),
                              metric = "n",
                              count_threshold = 40) {
  # Sum up unknown counts
  unknown_counts = unknown_counts[, .(n = sum(n)), by = id_cols]
  unknown_counts[, expected_read := FALSE]

  # Add additional metadata columns such as cell set to summed unknown counts
  meta_info = annotated_counts |> dplyr::distinct(across(all_of(c(id_cols, metadata_cols))))
  unknown_counts = dplyr::left_join(unknown_counts, meta_info, by = c(id_cols))

  # Calculate number of expected cell lines per well
  num_expected_lines = annotated_counts[expected_read == TRUE & is.na(cb_name)] |>
    dplyr::distinct(across(all_of(cell_line_cols))) |>
    dplyr::group_by(cell_set) |>
    dplyr::summarise(n_expected_lines = dplyr::n(), .groups = "drop")

  # Bind tables and calculate stats
  plate_well = dplyr::bind_rows(annotated_counts, unknown_counts) |>
    dplyr::left_join(num_expected_lines, by = "cell_set") |>
    dplyr::group_by(across(all_of(c(id_cols, metadata_cols)))) |>
    dplyr::summarise(
      # Total reads
      n_total_reads = sum(.data[[metric]], na.rm = TRUE),
      # Reads mapping to cell lines
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      # Fraction of reads mapping to cell lines
      fraction_expected_reads = n_expected_reads / n_total_reads,
      # Reads mapping to control barcodes
      n_cb_reads = sum(.data[[metric]][!is.na(cb_name)], na.rm = TRUE),
      # Median reads for control barcodes
      median_cb_reads = median(.data[[metric]][!is.na(cb_name)], na.rm = TRUE),
      # Fraction of reads mapping to control barcodes
      fraction_cb_reads = n_cb_reads / n_total_reads,
      # Number of cell lines with coverage above 40 reads
      n_lines_recovered = sum(.data[[metric]] >= count_threshold & is.na(cb_name) & expected_read == TRUE,
                              na.rm = TRUE),
      # Number of expected lines based on metadata
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      # Fraction of cell lines with coverage above count threshold
      fraction_cl_recovered = n_lines_recovered / n_expected_lines,
      # Ratio of control barcode reads to cell line reads
      cb_cl_ratio_well = n_cb_reads / n_expected_reads,
      .groups = "drop"
    )

  return(plate_well)
}

#' Calculate CB metrics
#'
#' This function calculates control metrics.
#'
#' @param normalized_counts A data frame containing normalized read counts and associated metrics.
#' @param cb_meta A data frame containing control barcode metadata.
#' @param group_cols A character vector specifying the grouping columns
#'                   (default: `c("pcr_plate", "pcr_well")`).
#' @param pseudocount An integer specifying the pseudocount to add to all read counts. Default: 20.
#' @return A data frame containing unique combinations of `group_cols` and the calculated metrics.
calculate_cb_metrics = function(normalized_counts, cb_meta,
                                group_cols = c("pcr_plate", "pcr_well"), pseudocount = 20) {
  # Filter cb_meta by dropping any control barcodes not used for normalization
  if ("cb_type" %in% colnames(cb_meta)) {
    dropped_cbs = cb_meta[cb_type != "well_norm"]

    if (nrow(dropped_cbs) > 0) {
      message("The following CBs in CB_meta are not used for normalization.")
      print(dropped_cbs)
      cb_meta = cb_meta[cb_type == "well_norm"]
    }
  }

  # Pull out control barcodes used for normalization and calculate some stats
  fit_stats = normalized_counts |>
    dplyr::semi_join(cb_meta, by = c("cb_ladder", "cb_name")) |>
    dplyr::group_by(across(all_of(group_cols)), across(any_of(c("cb_intercept", "log2_pseudovalue")))) |>
    dplyr::mutate(log2_norm_n = log2(n + pseudocount) + cb_intercept,
                  residual2 = (cb_log2_dose - log2_norm_n)^2,
                  squares2 = (cb_log2_dose - mean(cb_log2_dose))^2) |>
    dplyr::summarise(cb_mae = median(abs(cb_log2_dose - log2_norm_n)),
                     cb_r2 = 1 - sum(residual2) / sum(squares2),
                     cb_spearman = cor(cb_log2_dose, log2(n + pseudocount),
                                       method = "spearman", use = "pairwise.complete.obs"),
                     .groups = "drop")

  return(fit_stats)
}


# median cb cl ration plate and pert type - to do:fix thresholds
#' Calculate CB CL ratio of a plate
#'
#' This function calculates the cb cl ratio of a PCR plate. This takes in a read_stats df,
#' a cb_metrics df, and various thresholds.
#'
#' @param read_stats A dataframe with pcr_well level stats.
#' @param cb_metrics A dataframe with control barcode stats.
#' @param expected_reads_threshold A float
#' @param cb_threshold A integer indicating low control barcode reads.
#' @param cb_spearman_threshold A float indicating the Spearman correlation threshold.
#' @param cb_mae_threshold A integer indicating MAE threshold.
#' @return A data frame with CB CL ratio.
compute_cb_cl_ratio_plate = function(read_stats, cb_metrics,
                                     expected_reads_threshold = 0.8,
                                     cb_threshold = 40,
                                     cb_spearman_threshold = 0.8,
                                     cb_mae_threshold = 1) {
  cb_cl_ratio_plate = read_stats |>
    dplyr::left_join(cb_metrics, by = c("pcr_plate", "pcr_well"), suffix = c("", ".y")) |>
    dplyr::filter(median_cb_reads > cb_threshold,
                  fraction_expected_reads > expected_reads_threshold,
                  cb_spearman > cb_spearman_threshold & cb_mae < cb_mae_threshold) |>
    dplyr::group_by(pcr_plate, pert_type) |>
    dplyr::summarise(cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE), .groups = "drop")
  return(cb_cl_ratio_plate)
}

#' Compute Skew
#'
#' This function computes the skew, which measures the cumulative fraction of barcode read counts taken up by each cell
#' line. It is computed as the auc of that CDF function and has a range of (0.5,1).
#' A lower skew indicates a more even distribution of reads across cell lines.
#'
#' @param df A data frame containing the data to compute skew, generally annotated_counts.
#' @param group_cols A character vector specifying the column names to group by (default: `c("pcr_plate", "pcr_well")`).
#' @param metric A string indicating the column name of the metric to use for calculations (default: `"n"`).
#'
#' @return A data frame with one row per group and a column `skew` containing the computed skew (auc).
compute_skew = function(df, group_cols = c("pcr_plate", "pcr_well"), metric = "n") {
  result = df |>
    dplyr::group_by(across(all_of(group_cols))) |>
    arrange(desc(.data[[metric]])) |> # Sort by metric in descending order
    mutate(
      rank_fraction = row_number() / n(), # Calculate rank fraction of each cell line
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE)
      # Calculate cumulative fraction of reads
    ) |>
    summarise(
      skew = if (n() > 1) {
        # Compute auc
        sum((rank_fraction[-1] - rank_fraction[-n()]) *
          (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE)
      }
    ) |>
    dplyr::ungroup()
  return(result)
}


#' Generate ID column QC metrics table
#'
#' This function computes and combines various quality control (QC) metrics, such as read statistics, skew,
#' and control barcode metrics, grouped by a specified list of id_cols. For standard sequencing runs, this should always
#' be "pcr_plate" and "pcr_well".
#'
#' @param unknown_counts A data.table containing reads that did not map to any PRISM barcodes.
#' @param annotated_counts A data frame containing annotated read counts with metadata.
#' @param normalized_counts A data frame containing normalized read counts and associated metrics/metadata.
#' @param cb_meta
#' @param id_cols A character vector specifying the columns to group by for QC metric computations.
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line.
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered
#'                        (default: `40`).
#' @param cb_threshold A integer indicating low control barcode reads.
#' @param pseudocount An integer specifying the pseudocount to add to all read counts. Default: 20.
#' @param metadata_cols A character vector of column names of additional coluums to include.
#'   Default: 'c("pert_type", "pert_plate", "cell_set")'.
#' @return A data frame (`id_cols_table`) that combines QC metrics, including read statistics, skew, and
#'         control barcode metrics, grouped by the specified id_cols.
#'
#' @import dplyr
generate_id_cols_table = function(unknown_counts, annotated_counts, normalized_counts, cb_meta,
                                  id_cols, cell_line_cols,
                                  metadata_cols = c("pert_plate", "pert_type", "cell_set"),
                                  metric = "n",
                                  pseudocount = 20,
                                  count_threshold = 40,
                                  expected_reads_threshold = 0.8,
                                  cb_threshold = 40,
                                  cb_spearman_threshold = 0.8,
                                  cb_mae_threshold = 1) {
  # Columns to group read counts together to calculate stats - this should identify each PCR well
  message("Computing id_cols QC metric grouping over ",
          paste(id_cols, metadata_cols, collapse = ", "), "...")

  message("Calculating read count statistics...")
  read_stats = compute_read_stats(unknown_counts = unknown_counts,
                                  annotated_counts = annotated_counts,
                                  id_cols = id_cols,
                                  cell_line_cols = cell_line_cols,
                                  metadata_cols = metadata_cols,
                                  metric = metric,
                                  count_threshold = count_threshold)

  message("Calculating control barcode metrics...")
  cb_metrics = calculate_cb_metrics(normalized_counts, cb_meta,
                                    group_cols = id_cols,
                                    pseudocount = pseudocount)

  message("Calculating cb_cl_ratio for each PCR plate...")
  cb_cl_ratio_plate = compute_cb_cl_ratio_plate(read_stats, cb_metrics,
                                                expected_reads_threshold = expected_reads_threshold,
                                                cb_threshold = cb_threshold,
                                                cb_spearman_threshold = cb_spearman_threshold,
                                                cb_mae_threshold = cb_mae_threshold)

  message("Calculating read count skew...")
  skew = compute_skew(annotated_counts, group_cols = id_cols, metric = metric)

  message("Assembling id_cols_table...")
  id_cols_table = read_stats |>
    dplyr::left_join(cb_cl_ratio_plate, by = c("pcr_plate", "pert_type")) |>
    dplyr::left_join(skew, by = id_cols) |>
    dplyr::left_join(cb_metrics, by = id_cols)

  return(id_cols_table)
}

#' Generate QC Flags for ID Columns in Plate Cell Data
#'
#' This function processes a data frame containing plate cell measurements and assigns quality control (QC) flags
#' based on several thresholds. The flags are determined in a sequential order (i.e., only the first applicable
#' flag per well is recorded). The function returns a data frame with a unique record of flagged wells, identified
#' by the specified grouping columns.
#'
#' @param id_cols_table A data frame containing plate cell measurements. It must include the following columns:
#'   \code{median_cb_reads}, \code{fraction_expected_reads}, \code{cb_mae}, \code{cb_spearman},
#'   \code{cb_cl_ratio_well}, and \code{pert_type}.
#' @param group_cols A character vector specifying the columns used to uniquely identify each well. Default is
#'   \code{c("pcr_plate", "pcr_well", "pert_type")}.
#' @param contamination_threshold Numeric threshold for contamination based on the fraction of expected reads.
#'   Default is 0.8.
#' @param cb_mae_threshold Numeric threshold for the CB mean absolute error (MAE). Default is 1.
#' @param cb_spearman_threshold Numeric threshold for the CB Spearman correlation. Default is 0.8.
#' @param cb_cl_ratio_low_negcon Numeric threshold for the lower bound of the CB/CL ratio in negative control wells.
#'   Default is 0.
#' @param cb_cl_ratio_high_negcon Numeric threshold for the upper bound of the CB/CL ratio in negative control wells.
#'   Default is 2.
#' @param cb_cl_ratio_low_poscon Numeric threshold for the lower bound of the CB/CL ratio in positive control wells.
#'   Default is 0.5.
#' @param cb_cl_ratio_high_poscon Numeric threshold for the upper bound of the CB/CL ratio in positive control wells.
#'   Default is 2.
#' @param well_reads_threshold Numeric threshold for well reads; wells with \code{median_cb_reads} below the logarithm
#'   of this value are flagged as "well_reads". Default is 100.
#'
#' @return A data frame containing the unique flagged wells. The returned data frame includes the columns specified
#'   in \code{group_cols} and a \code{qc_flag} column that indicates the first QC flag applied to each well.
#'
#' @import dplyr
id_cols_qc_flags <- function(id_cols_table,
                             group_cols = c("pcr_plate", "pcr_well", "pert_type", "pert_plate", "cell_set"),
                             contamination_threshold = 0.8,
                             cb_mae_threshold = 1,
                             cb_spearman_threshold = 0.8,
                             cb_cl_ratio_low_negcon = 0,
                             cb_cl_ratio_high_negcon = 2,
                             cb_cl_ratio_low_poscon = 0.5,
                             cb_cl_ratio_high_poscon = 2,
                             well_reads_threshold = 40) {
  # Add a qc_flag column using case_when (conditions are checked in order)
  qc_table <- id_cols_table %>%
    mutate(qc_flag = case_when(
      median_cb_reads < well_reads_threshold ~ "well_reads",
      fraction_expected_reads < contamination_threshold ~ "contamination",
      cb_mae > cb_mae_threshold | cb_spearman < cb_spearman_threshold ~ "cb_linearity",
      pert_type == "trt_poscon" & (cb_cl_ratio_well < cb_cl_ratio_low_poscon |
        cb_cl_ratio_well > cb_cl_ratio_high_poscon) ~ "cb_cl_ratio",
      pert_type == "ctl_vehicle" & (cb_cl_ratio_well < cb_cl_ratio_low_negcon |
        cb_cl_ratio_well > cb_cl_ratio_high_negcon) ~ "cb_cl_ratio",
      TRUE ~ NA_character_
    ))

  # Extract flagged wells (keeping only the first flag per well)
  flagged_all <- qc_table %>%
    filter(!is.na(qc_flag)) %>%
    select(all_of(group_cols), qc_flag) %>%
    unique()

  # Return both outputs
  return(flagged_all)
}

#' Generate Pool Well QC Table
#'
#' This function flags pool/well combinations in control wells based on the variability of cell line
#' measurements relative to the pool/well median. It calculates the median normalized count for each
#' pool/well group, computes the absolute difference between each cell line's normalized count and the
#' group median, and then determines the fraction of outliers. Wells with a fraction of outliers exceeding
#' the specified threshold are flagged as having "pool_well_outliers".
#'
#' @param normalized_counts A data frame containing normalized count data. Required columns include:
#'   \code{cb_name}, \code{pcr_plate}, \code{pcr_well}, \code{pert_type}, \code{pool_id},
#'   \code{log2_normalized_n}, \code{lua}, \code{depmap_id}, and \code{cell_set}.
#' @param pool_well_delta_threshold A numeric threshold specifying the minimum absolute difference from
#'   the pool/well median for a cell line to be considered an outlier (default: 5).
#' @param pool_well_fraction_threshold A numeric threshold specifying the minimum fraction of outlier cell
#'   lines required for the well to be flagged (default: 0.4).
#'
#' @return A data frame with an added \code{qc_flag} column that indicates pool/well combinations flagged
#'   as having excessive outliers.
#'
#' @import dplyr
generate_pool_well_qc_table <- function(normalized_counts,
                                        id_cols = c("pcr_plate", "pcr_well"),
                                        pool_well_delta_threshold = 5,
                                        pool_well_fraction_threshold = 0.4) {
  ### POOL_WELL_OUTLIERS
  ## Flag pool/well combinations based on the fraction of cell lines in a pool + well that are some distance
  # from the pool + well median.
  ## This is in control wells only
  # To do: Comments mention negcons!!! But this is doing all wells?!?!?!
  pool_well_outliers = normalized_counts |>
    # Get the median value of the pool in each well
    dplyr::group_by(across(all_of(c(id_cols, "pert_type", "pool_id")))) |>
    dplyr::mutate(pool_well_median = median(log2_normalized_n, na.rm = TRUE),
                  delta_from_pool_well_median = abs(log2_normalized_n - pool_well_median)) |>
    dplyr::summarise(n_outliers = sum(delta_from_pool_well_median > pool_well_delta_threshold, na.rm = TRUE),
                     n_cell_lines = dplyr::n(),
                     .groups = "drop") |>
    dplyr::mutate(fraction_outliers = n_outliers / n_cell_lines,
                  qc_flag = if_else(fraction_outliers > pool_well_fraction_threshold,
                                    "pool_well_outliers", NA_character_))
}

#' Compute variance decomposition
#'
#' This function calculates the variabilty of cell lines within plate and
#' the variability of the pools.
#'
#' @param normalized_counts
#' @param metric
#' @param negcon
#' @param cell_line_cols
#' @param id_cols
#' @return A data frame
compute_variance_decomposition <- function(normalized_counts, metric = "n", negcon = "ctl_vehicle",
                                           cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
                                           id_cols = c("pcr_plate", "pcr_well")) {
  # Add "pert_plate" to id_cols if not already present
  id_cols <- c(id_cols, "pert_plate")

  # Add pool_id annotations to control pools
  df <- normalized_counts %>%
    dplyr::mutate(pool_id=ifelse(!is.na(cb_name), "CTLBC", pool_id))

  # Compute variance of cell line fractions
  var_log_fline <- df %>%
    filter(pert_type==negcon) %>%
    dplyr::group_by(across(c(all_of(id_cols), cell_set))) %>%
    dplyr::mutate(tot_counts=sum(.data[[metric]]+1)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(across(c(all_of(id_cols), all_of(cell_line_cols), cb_name))) %>%
    dplyr::mutate(fcell_line=sum(.data[[metric]]+1)/tot_counts) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(across(c(all_of(cell_line_cols), cb_name, pcr_plate, pert_plate))) %>%
    dplyr::summarise(var_log2_fline=var(log2(fcell_line)),
                      median_log2_fline=median(log2(fcell_line)),
                      mad_log2_fline=mad(log2(fcell_line)),
                      mean_log2_fline=mean(log2(fcell_line))) %>%
    dplyr::ungroup()


  # Compute variance of cell line fractions in pools
  var_log_cl_in_pool <- df %>%
      filter(pert_type==negcon) %>%
      dplyr::group_by(across(c(all_of(id_cols), cell_set, pool_id))) %>%
      dplyr::mutate(tot_counts=sum(.data[[metric]]+1)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(across(c(all_of(id_cols), all_of(cell_line_cols), cb_name))) %>%
      dplyr::mutate(fcl_in_pool=sum(.data[[metric]]+1)/tot_counts) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(across(c(all_of(cell_line_cols), cb_name, pcr_plate, pert_plate))) %>%
      dplyr::summarise(var_log2_fcl_in_pool=var(log2(fcl_in_pool)),
                        mad_log2_fcl_in_pool=mad(log2(fcl_in_pool)),
                        median_log2_fcl_in_pool=median(log2(fcl_in_pool)),
                        mean_log2_fcl_in_pool=mean(log2(fcl_in_pool))) %>%
      dplyr::ungroup()

  # Compute fraction of reads in pools
  pwise_negcon_stats <- df %>%
      filter(pert_type==negcon) %>%
      dplyr::group_by(cell_set, across(all_of(id_cols))) %>%
      dplyr::mutate(tot_counts=sum(.data[[metric]]+1)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(across(c(all_of(id_cols), cell_set, pool_id))) %>%
      dplyr::summarise(
          frac_reads=sum(.data[[metric]]+1)/tot_counts) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()


  # Compute variance of fraction of reads in pools
  var_log_fpool <- pwise_negcon_stats %>%
      dplyr::group_by(cell_set, pcr_plate, pert_plate,
                      pool_id) %>%
      dplyr::summarise(var_log2_frac_pool_reads=var(log2(frac_reads)),
                        mad_log2_frac_pool_reads=mad(log2(frac_reads)),
                        median_log2_frac_pool_reads=median(log2(frac_reads)),
                        mean_log2_frac_pool_reads=mean(log2(frac_reads))) %>%
      dplyr::ungroup()


  # Join all variance components together
  var_decomp <- dplyr::left_join(var_log_cl_in_pool, var_log_fpool,
                                      by=c("pool_id", "cell_set", "pcr_plate", "pert_plate")) %>%
      dplyr::left_join(var_log_fline) %>%
    dplyr::group_by(cell_set, pool_id, pcr_plate, pert_plate) %>%
    dplyr::summarise(across(where(is.numeric), function(x) median(x, na.rm = TRUE)))

  return(var_decomp)
}

#' Create contamination tables
#'
#' This creates three tables investigating contaminants.
#'
#' @param prism_barcode_counts
#' @param unknown_barcode_counts
#' @param cell_set_and_pool_meta
#' @param cell_line_meta
#' @param cb_meta
#' @param sample_meta
#' @return A list of data frames
compute_contamination_qc_tables <- function(prism_barcode_counts,
                                            unknown_barcode_counts,
                                            cell_set_and_pool_meta,
                                            cell_line_meta,
                                            cb_meta,
                                            sample_meta,
                                            id_cols = c("pcr_plate", "pcr_well")) {
  # --- 0. Ensure all required columns are present ---
  prism_barcode_counts_all <- prism_barcode_counts %>%
    dplyr::left_join(sample_meta %>% select(pcr_plate, pcr_well, pert_plate), by = id_cols)

  unknown_barcode_counts_all <- unknown_barcode_counts %>%
    dplyr::left_join(sample_meta %>% select(pcr_plate, pcr_well, pert_plate), by = id_cols)

  # --- 1. Create total counts df with known and unknown barcodes ---
  total_counts <- bind_rows(prism_barcode_counts_all, unknown_barcode_counts_all)

  # --- 2. Compute total counts for each well ---
  total_counts_by_well <- total_counts %>%
    group_by(pcr_plate, pcr_well, pert_plate) %>%
    summarise(well_count = sum(n), .groups = "drop")

  # --- 3. Get a list of expected reads ---
  expected_cell_lines <- cell_set_and_pool_meta %>%
    left_join(cell_line_meta, by = c("depmap_id", "lua")) %>%
    pull(forward_read_barcode) %>%
    unique()

  expected_controls <- cb_meta %>%
    pull(forward_read_barcode) %>%
    unique()

  unexpected_cell_lines <- cell_line_meta %>%
    anti_join(cell_set_and_pool_meta, by = "depmap_id") %>%
    pull(forward_read_barcode) %>%
    unique()

  # --- 4. Annotate the counts with the read type ---
  total_counts_with_read_type <- total_counts %>%
    mutate(
      read_type = case_when(
        forward_read_barcode %in% expected_cell_lines ~ "expected_cell_line",
        forward_read_barcode %in% expected_controls ~ "control",
        forward_read_barcode %in% unexpected_cell_lines ~ "unexpected_cell_line",
        forward_read_barcode == "unknown_low_abundance_barcode" ~ "unknown_low_abundance",
        TRUE ~ "unknown_barcode" # otherwise() is equivalent to TRUE ~
      )
    )

  # --- 5. Annotate with sample metadata ---
  total_counts_with_read_type_annotated <- total_counts_with_read_type %>%
    left_join(sample_meta, by = c("pcr_plate", "pcr_well", "pert_plate"))

  # --- 6. Join with total well counts and aggregate ---
  well_counts_by_read_type <- total_counts_with_read_type_annotated %>%
    left_join(total_counts_by_well, by = c("pcr_plate", "pcr_well", "pert_plate")) %>%
    group_by(pcr_plate, pcr_well, pert_plate, read_type) %>%
    summarise(
      n = sum(n),
      well_count = first(well_count), # pl.first() is equivalent to first()
      .groups = "drop"
    )

  # --- 7. Compute the fraction of each read type per well ---
  fraction_read_type_by_well <- well_counts_by_read_type %>%
    mutate(fraction_well_reads = n / well_count)

  # --- 8. Compute the mean of the fractions per plate ---
  fraction_read_type_by_plate <- fraction_read_type_by_well %>%
    group_by(pcr_plate, pert_plate, read_type) %>%
    summarise(mean_fraction_reads = mean(fraction_well_reads), .groups = "drop")

  # --- 9. Filter the total_counts_with_read_type_annotated to include only relevant columns ---
  total_counts_by_read_type <- total_counts_with_read_type_annotated %>%
    select(c("pcr_plate","pert_plate","pcr_well","forward_read_barcode","read_type","n"))

  return(list(fraction_read_type_by_well = fraction_read_type_by_well,
              fraction_read_type_by_plate = fraction_read_type_by_plate,
              total_counts_by_read_type = total_counts_by_read_type))
}
