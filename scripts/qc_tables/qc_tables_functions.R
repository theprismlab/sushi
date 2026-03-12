options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(PRROC)
library(dplyr)


# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

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
#'
#' @import dplyr
compute_skew <- function(df, group_cols = c("pcr_plate", "pcr_well", "pert_plate"), metric = "n") {
  result <- df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    arrange(desc(.data[[metric]])) %>% # Sort by metric in descending order
    mutate(
      rank_fraction = row_number() / n(), # Calculate rank fraction of each cell line
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE) # Calculate cumulative fraction of reads
    ) %>%
    summarise(
      skew = if (n() > 1) {
        # Compute auc
        sum((rank_fraction[-1] - rank_fraction[-n()]) *
          (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE)
      }
    ) %>%
    dplyr::ungroup()
  return(result)
}

#' Compute expected cell lines
#'
#' This function calculates the number of unique expected cell lines for each cell_set.
#'
#' @param cell_set_meta A data frame containing metadata for cell sets.
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line.
#'
#' @return A data frame with one row per cell_set and a column `n_expected_lines` indicating the number of unique expected cell lines in each set.
#'
#' @import dplyr
compute_expected_lines <- function(cell_set_meta, cell_line_cols) {
  cell_set_meta %>%
    dplyr::distinct(cell_set, across(all_of(cell_line_cols))) %>%
    dplyr::group_by(cell_set) %>%
    dplyr::summarise(n_expected_lines = n(), .groups = "drop")
}

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
    dplyr::left_join(cb_metrics, by = c("pcr_plate", "pcr_well"), suffix = c("", ".y")) %>%
    dplyr::filter(median_cb_reads > cb_threshold,
                  fraction_expected_reads > expected_reads_threshold,
                  cb_spearman > cb_spearman_threshold & cb_mae < cb_mae_threshold) %>%
    dplyr::group_by(pcr_plate, pert_type) %>%
    dplyr::summarise(cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE), .groups = "drop")
  return(cb_cl_ratio_plate)
}

# TABLE GENERATION FUNCTION ----------


# CELL LINE BY PLATE (PCR_PLATE, CELL LINE) ----------

#' Compute error rate
#'
#' This function calculates the error rate using receiver operating characteristic (ROC) curve data.
#' It quantifies the overlap or misclassification between positive controls and negative controls
#' within grouped data, providing a measure of assay quality. A high error rate indicates large overlap.
#'
#' @param df A data frame containing data for analysis, including metrics and annotations for negative and positive controls.
#' @param metric A string specifying the column name of the metric to use for computing error rate (default: `"log2_normalized_n"`).
#' @param group_cols A character vector specifying the grouping columns (default: `c("depmap_id", "pcr_plate")`).
#' @param negcon A string specifying the `pert_type` value that identifies negative control samples (default: `"ctl_vehicle"`).
#'   These represent baseline or no-treatment conditions and are expected to have lower values.
#' @param poscon A string specifying the `pert_type` value that identifies positive control samples (default: `"trt_poscon"`).
#'   These represent treated conditions and are expected to have higher values.
#'
#' @return A data frame with one row per group and a column `error_rate` containing the computed error rate.
#'
#' @details
#' - **How the Error Rate is Calculated**:
#'   - The function filters the data to include only rows where `pert_type` matches `negcon` or `poscon`.
#'   - For each group defined by `group_cols`, the function computes an ROC curve using the specified `metric`.
#'   - The `metric` values are treated as scores, and the `pert_type` column is used to assign class labels:
#'     - `negcon` samples are treated as the negative class (class 0).
#'     - `poscon` samples are treated as the positive class (class 1).
#'   - The error rate is calculated as the minimum misclassification rate, which is derived from the ROC curve:
#'     - \(\text{Error Rate} = \min\left(\text{True Positive Rate} + 1 - \text{True Negative Rate}\right) / 2\).
#'
#' @import dplyr
#' @import PRROC
compute_error_rate <- function(df, metric = "log2_normalized_n", group_cols = c("depmap_id", "pcr_plate", "pert_plate"),
                               negcon = "ctl_vehicle", poscon = "trt_poscon", contains_poscon = TRUE) {

  if (contains_poscon) {
    message("Computing error rate using ", negcon, " and ", poscon, ".....")
    message("Grouping by ", paste(group_cols, collapse = ", "), ".....")
    result <- df %>%
      dplyr::filter(
        pert_type %in% c(negcon, poscon),
        is.finite(.data[[metric]]),
        !is.na(pool_id)
      ) %>%
      dplyr::group_by(across(all_of(group_cols))) %>%
      dplyr::summarise(
        error_rate = {
          roc_data <- PRROC::roc.curve(
            scores.class0 = .data[[metric]],
            weights.class0 = pert_type == negcon,
            curve = TRUE
          )
          min(roc_data$curve[, 1] + 1 - roc_data$curve[, 2]) / 2
        },
      ) %>%
      dplyr::ungroup()
    return(result)
  } else {
    print("No positive controls found. Unable to calculate error rate.")
  }
}

#' Compute control median and MAD
#'
#' This function calculates the median and median absolute deviation (MAD) for negative controls
#' and positive controls in both raw and normalized data. It also computes false sensitivity probabilities
#' based on the MAD values for the negative controls.
#'
#' @param df A data frame containing normalized and raw data along with annotations for control types.
#' @param group_cols A character vector specifying the columns to group by, along with `pert_type` (default: `c("depmap_id", "pcr_plate")`).
#' @param negcon A string specifying the `pert_type` value for negative controls (default: `"ctl_vehicle"`).
#' @param poscon A string specifying the `pert_type` value for positive controls (default: `"trt_poscon"`).
#'
#' @return A data frame containing medians and MADs for raw and normalized data for both control types.
#' Additional columns include false sensitivity probabilities at thresholds -1 (`false_sensitivity_probability_50`)
#' and -2 (`false_sensitivity_probability_25`).
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom stats mad pnorm
compute_ctl_medians_and_mad <- function(df, group_cols = c("depmap_id", "lua", "pcr_plate", "pert_plate"),
                                        negcon = "ctl_vehicle", poscon = "trt_poscon", pseudocount = 20) {
  message("Adding median and MAD values for ", negcon, " and ", poscon, " if it exists.....")
  message("Computing falses sensitivity probability for ", negcon, ".....")
  # Group and compute medians/MADs
  result <- df %>%
    dplyr::filter(pert_type %in% c(negcon, poscon)) %>%
    dplyr::group_by(across(all_of(c(group_cols, "pert_type", "day")))) %>%
    dplyr::summarise(
      median_log_normalized = median(log2_normalized_n),
      n_replicates = n(),
      mad_log_normalized = mad(log2_normalized_n),
      median_raw = median(n),
      mad_raw = mad(log2(n + pseudocount))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      names_from = pert_type,
      values_from = c(median_log_normalized, mad_log_normalized, median_raw, mad_raw, n_replicates),
      names_sep = "_"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      false_sensitivity_probability_50 = pnorm(
        -1,
        sd = .data[[paste0("mad_log_normalized_", negcon)]] * sqrt((1 / .data[[paste0("n_replicates_", negcon)]] * pi / 2) + 1)
      ),
      false_sensitivity_probability_25 = pnorm(
        -2,
        sd = .data[[paste0("mad_log_normalized_", negcon)]] * sqrt((1 / .data[[paste0("n_replicates_", negcon)]] * pi / 2) + 1)
      )
    )
  return(result)
}

#' Compute log fold change of positive controls
#'
#' This function calculates the log fold change (LFC) between positive controls (`poscon`) and negative controls (`negcon`)
#' based on median values from normalized and raw data.
#'
#' @param df A data frame containing median values for both normalized and raw data, with separate columns for positive and negative controls.
#' @param negcon A string specifying the prefix for columns corresponding to the negative control (default: `"ctl_vehicle"`).
#' @param poscon A string specifying the prefix for columns corresponding to the positive control (default: `"trt_poscon"`).
#'
#' @return A data frame with additional columns:
#' - `lfc_normalized`: Log fold change for normalized data.
#' - `lfc_raw`: Log fold change for raw data.
#'
#' @import dplyr
compute_control_lfc <- function(df, negcon = "ctl_vehicle", poscon = "trt_poscon",
                                grouping_cols = c("depmap_id", "pcr_plate", "pert_plate"), 
                                contains_poscon = TRUE) {
  message("Computing log fold change for ", negcon, " and ", poscon, ".....")
  if (contains_poscon) {
    result <- df %>%
      dplyr::mutate(
        lfc_trt_poscon = .data[[paste0("median_log_normalized_", poscon)]] -
          .data[[paste0("median_log_normalized_", negcon)]],
        lfc_raw_trt_poscon = .data[[paste0("median_raw_", poscon)]] -
          .data[[paste0("median_raw_", negcon)]]
      ) %>%
      dplyr::select(all_of(grouping_cols), lfc_trt_poscon, lfc_raw_trt_poscon) %>%
      dplyr::mutate(viability_trt_poscon = 2^lfc_trt_poscon)
    return(result)
  } else {
    print("No positive controls found. Unable to calculate log fold change.")
  }
}

#' Compute cell line fractions
#'
#' This function computes the total reads and fraction of reads contributed by each cell line within groups defined by specified columns.
#'
#' @param df A data frame containing read data for cell lines, including a metric for read counts.
#' @param metric A string specifying the column name of the metric to use for calculations (default: `"n"`).
#' @param grouping_cols A character vector specifying the columns to group by (default: `c("pcr_plate", "depmap_id")`).
#'
#' @return A data frame with the following additional columns for each group:
#' - `total_reads`: The total reads per group.
#' - `fraction_of_reads`: The fraction of total reads contributed by each entry.
#'
#' @import dplyr
compute_cl_fractions <- function(df, metric = "n", grouping_cols = c("pcr_plate", "depmap_id", "pert_plate")) {
  print(paste0("Computing cell line fractions for ", metric, "....."))
  result <- df %>%
    dplyr::group_by(across(all_of(grouping_cols))) %>%
    dplyr::summarise(
      total_reads = sum(.data[[metric]], na.rm = TRUE), # Total reads per group
      fraction_of_reads = sum(.data[[metric]], na.rm = TRUE) / sum(.data[[metric]], na.rm = TRUE) # Fraction of reads for each entry
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(grouping_cols), total_reads, fraction_of_reads)
  return(result)
}

#' Compute median number of biological replicates in treatments
#'
#' Actions:
#' Grab normalized counts, filter for treatments
#' Group by cell lines + sig cols + plate and count number of bio reps,
#' Group by cell lines + plate and get median number of bio reps
#'
#' @param norm_counts A dataframe of filtered nromalized counts.
#' @param cell_line_cols A vector of columns describing cell lines.
#' @param sig_cols A vector of columns describing treatment profiles.
#' @return A dataframe.
compute_med_trt_bio_rep = function(norm_counts, cell_line_cols, sig_cols) {
  med_trt_bio_reps = norm_counts |>
    dplyr::filter(pert_type == "trt_cp") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(unique(c(cell_line_cols, sig_cols, "pert_plate"))))) |>
    dplyr::summarise(num_trt_bio_reps = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(unique(c(cell_line_cols, "pert_plate"))))) |>
    dplyr::summarise(med_num_trt_bio_reps = median(num_trt_bio_reps), .groups = "drop")

  return(med_trt_bio_reps)
}

#' Generate cell plate table
#'
#' This function generates a comprehensive QC table for cell_lines +pcr plates by computing and merging various QC metrics, including medians, MADs, error rates, log fold changes (LFC), and cell line fractions.
#'
#' @param normalized_counts A data frame containing normalized read counts and associated metadata.
#' @param filtered_counts A data frame containing filtered read counts for the same dataset.
#' @param cell_line_cols A string of comma-separated column names that define a unique cell line.
#'
#' @return A data frame (`plate_cell_table`) that merges various QC metrics grouped by cell lines and plates, including:
#' - Control medians and MADs for normalized and raw data.
#' - Error rates based on positive and negative control separation.
#' - Log fold changes (LFC) for positive controls.
#' - Fractions of reads contributed by each cell line.
#'
#' @import dplyr
generate_plate_cell_table = function(normalized_counts, sample_meta,
                                     cell_line_cols, sig_cols,
                                     pseudocount = 20, contains_poscon = TRUE,
                                     poscon = "trt_poscon", negcon = "ctl_vehicle",
                                     nc_variability_threshold = 1, error_rate_threshold = 0.05,
                                     pc_viability_threshold = 0.25, nc_raw_count_threshold = 40) {

  # QC outputs require more columns in addition to cell_line_cols.
  cell_line_plate_grouping = unique(c(cell_line_cols, "pcr_plate", "pert_plate", "project_code", "day"))
  message("Computing cell + plate QC metrics grouping by ", paste(cell_line_plate_grouping, collapse = ", "), "...")

  # Check that all columns are present in normalized_counts
  missing_cols = setdiff(cell_line_plate_grouping, colnames(normalized_counts))
  if (length(missing_cols)) {
    stop("The following columns are missing from normalized counts: ", paste(missing_cols, collapse = ", "))
  }

  # Compute control medians and MAD
  message("Calculating medians and MAD for controls ...")
  medians_and_mad = compute_ctl_medians_and_mad(
    df = normalized_counts,
    group_cols = cell_line_plate_grouping,
    negcon = negcon,
    poscon = poscon,
    pseudocount = pseudocount
  )

  # Decide: Compute number of expected control wells from the sample meta
  message("Determining the expected number of controls on each PCR plate ...")
  n_expected_controls = sample_meta |>
    dplyr::filter(pert_type %in% c(negcon, poscon)) |>
    dplyr::group_by(pert_plate, pcr_plate, pert_type) |>
    dplyr::summarize(unique_bio_rep = n_distinct(bio_rep), .groups = "drop") |>
    tidyr::pivot_wider(names_from = pert_type, values_from = unique_bio_rep, names_prefix = "n_expected_")

  # Decide: Compute cell line fractions per plate
  message("Computing cell line fractions ...")
  cell_line_fractions = compute_cl_fractions(
    df = normalized_counts,
    grouping_cols = cell_line_plate_grouping
  )

  # Calc median number of bio reps across the treatments for each cell line + plate
  message("Computing median number of bio_reps across the treatments ...")
  med_trt_bio_reps = compute_med_trt_bio_rep(norm_counts = normalized_counts,
                                             cell_line_cols = cell_line_plate_grouping,
                                             sig_cols = sig_cols)

  # Join tables together
  message("Left joining tables ...")
  plate_cell_table = medians_and_mad |>
    dplyr::left_join(n_expected_controls, by = c("pcr_plate", "pert_plate")) |>
    dplyr::left_join(cell_line_fractions, by = cell_line_plate_grouping) |>
    dplyr::left_join(med_trt_bio_reps, by = cell_line_plate_grouping)

  if (contains_poscon) {
    # Compute error rate
    message("Detected poscons - computing error rates ...")
    error_rates = compute_error_rate(
      df = normalized_counts,
      metric = "log2_normalized_n",
      group_cols = cell_line_plate_grouping,
      negcon = negcon,
      poscon = poscon,
      contains_poscon = contains_poscon
    )

    # Compute poscon LFC
    message("Detected poscons - computing poscon l2fcs ...")
    poscon_lfc = compute_control_lfc(
      df = medians_and_mad,
      negcon = negcon,
      poscon = poscon,
      grouping_cols = cell_line_plate_grouping,
      contains_poscon = contains_poscon
    )

    # Left join poscon QC metrics, calculate fraction of expected ctrls,
    # and check qc values against thresholds
    message("Detected poscons - joining additional columns and checking against QC thresholds ...")
    plate_cell_table = plate_cell_table |>
      dplyr::left_join(error_rates, by = cell_line_plate_grouping) |>
      dplyr::left_join(poscon_lfc, by = cell_line_plate_grouping) |>
      dplyr::mutate(fraction_expected_negcon = .data[[paste0("n_replicates_", negcon)]] /
                      .data[[paste0("n_expected_", negcon)]],
                    fraction_expected_poscon = .data[[paste0("n_replicates_", poscon)]] /
                      .data[[paste0("n_expected_", poscon)]],
                    qc_pass = error_rate < error_rate_threshold &
                      .data[[paste0("viability_", poscon)]] < pc_viability_threshold &
                      .data[[paste0("median_raw_", negcon)]] > nc_raw_count_threshold &
                      .data[[paste0("mad_log_normalized_", negcon)]] < nc_variability_threshold)

  } else {
    # Calculate fraction of expected ctrls and check qc values against thresholds
    message("Checking against QC thresholds ...")
    plate_cell_table = plate_cell_table |>
      dplyr::mutate(fraction_expected_negcon = .data[[paste0("n_replicates_", negcon)]] /
                      .data[[paste0("n_expected_", negcon)]],
                    qc_pass = .data[[paste0("median_raw_", negcon)]] > nc_raw_count_threshold &
                      .data[[paste0("mad_log_normalized_", negcon)]] < nc_variability_threshold)
  }

  # Determine if cell line passes on a pert plate level
  message("Check qc pass rate at pert plate level ...")
  plate_cell_table = plate_cell_table |>
    dplyr::mutate(n_passing_med_num_trt_reps = ifelse(qc_pass, med_num_trt_bio_reps, 0)) |>
    dplyr::group_by(across(all_of(c(cell_line_cols, "pert_plate")))) |>
    dplyr::mutate(qc_pass_pert_plate = sum(n_passing_med_num_trt_reps) > 1) |>
    dplyr::ungroup()

  return(plate_cell_table)
}

# QC FLAG FUNCTIONS ----------

#' Compute QC Flags for Plate Cell Data
#'
#' This function processes a plate cell table by applying a series of quality control
#' thresholds to assign a flag to each well. The flag indicates the first QC criterion that
#' a well fails, based on the provided thresholds.
#'
#' @param plate_cell_table A data frame containing plate cell data. It must include the columns:
#'   \code{mad_log_normalized_ctl_vehicle}, \code{error_rate}, \code{viability_trt_poscon}, and
#'   \code{median_raw_ctl_vehicle}.
#' @param nc_variability_threshold A numeric threshold for negative control variability
#'   (default: 1).
#' @param error_rate_threshold A numeric threshold for the error rate (default: 0.05).
#' @param pc_viability_threshold A numeric threshold for positive control viability
#'   (default: 0.25).
#' @param nc_raw_count_threshold A numeric threshold for negative control raw counts, where the
#'   raw count is compared to the log of this value (default: 40).
#'
#' @return A data frame identical to \code{plate_cell_table} with an additional column \code{qc_flag}
#'   that indicates the first QC flag applicable for each well.
#'
#' @import dplyr

plate_cell_qc_flags <- function(plate_cell_table,
                                nc_variability_threshold = 1,
                                error_rate_threshold = 0.05,
                                pc_viability_threshold = 0.25,
                                nc_raw_count_threshold = 40,
                                contains_poscon = TRUE) {
  # Add a qc_flag column using case_when (conditions are checked in order)
  if (contains_poscon) {
    qc_table <- plate_cell_table |>
      mutate(qc_flag = case_when(
        mad_log_normalized_ctl_vehicle > nc_variability_threshold ~ "nc_variability",
        error_rate > error_rate_threshold ~ "error_rate",
        viability_trt_poscon > pc_viability_threshold ~ "pc_viability",
        median_raw_ctl_vehicle < log(nc_raw_count_threshold) ~ "nc_raw_count",
        TRUE ~ NA_character_
      ))
  } else {
    qc_table <- plate_cell_table %>%
      mutate(qc_flag = case_when(
        mad_log_normalized_ctl_vehicle > nc_variability_threshold ~ "nc_variability",
        median_raw_ctl_vehicle < log(nc_raw_count_threshold) ~ "nc_raw_count",
        TRUE ~ NA_character_
      ))
  }
  return(qc_table)
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
generate_pool_well_qc_table <- function(normalized_counts, pool_well_delta_threshold = 5, pool_well_fraction_threshold = 0.4) {
  ### POOL_WELL_OUTLIERS
  ## Flag pool/well combinations based on the fraction of cell lines in a pool + well that are some distance from the pool + well median.
  ## This is in control wells only
  # To do: this flags all wells not just the control wells
  pool_well_outliers = normalized_counts |>
    # Consider only cell line barcodes
    dplyr::filter(is.na(cb_name)) |>
    # Get the median value of the pool in each well
    dplyr::group_by(pcr_plate, pcr_well, pert_type, pool_id) |>
    dplyr::mutate(pool_well_median = median(log2_normalized_n, na.rm = TRUE),
                  n_cell_lines = dplyr::n_distinct(lua, depmap_id, cell_set),
                  delta_from_pool_well_median = abs(log2_normalized_n - pool_well_median)) |>
    dplyr::summarise(n_outliers = sum(delta_from_pool_well_median > pool_well_delta_threshold, na.rm = TRUE),
                     n_cell_lines = max(n_cell_lines),
                     .groups = "drop") |>
    dplyr::mutate(fraction_outliers = n_outliers / n_cell_lines,
                  qc_flag = if_else(fraction_outliers > pool_well_fraction_threshold,
                                    "pool_well_outliers", NA_character_))
}

# Get the qc parameters from the qc_params json file
load_thresholds_from_json <- function(json_file_path) {
  # Load required package
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.")
  }

  # Read JSON file into a list
  params <- jsonlite::fromJSON(json_file_path)

  # Convert each value to numeric
  numeric_params <- lapply(params, as.numeric)

  return(numeric_params)
}

# PCR PLATE FLAGS
generate_pcr_plate_qc_flags_table <- function(plate_cell_table, fraction_expected_controls, contains_poscon = TRUE) {
  # Add a qc_flag when either fraction_expected_poscon or fraction_expected_negcon is below the threshold
  if (contains_poscon) {
  table <- plate_cell_table %>%
    dplyr::select(fraction_expected_poscon, fraction_expected_negcon, pcr_plate, pert_plate) %>%
    unique() %>%
    dplyr::mutate(
      qc_flag = dplyr::case_when(
        fraction_expected_poscon < fraction_expected_controls ~ "fraction_expected_controls",
        fraction_expected_negcon < fraction_expected_controls ~ "fraction_expected_controls",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(qc_flag))
  } else {
    table <- plate_cell_table %>%
      dplyr::select(fraction_expected_negcon, pcr_plate, pert_plate) %>%
      unique() %>%
      dplyr::mutate(
        qc_flag = dplyr::case_when(
          fraction_expected_negcon < fraction_expected_controls ~ "fraction_expected_controls",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(qc_flag))
  }
  return(table)
}
