options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(PRROC)
library(dplyr)

#' Compute error rate
#'
#' This function calculates the error rate using receiver operating characteristic (ROC) curve data.
#' It quantifies the overlap or misclassification between positive controls and negative controls
#' within grouped data, providing a measure of assay quality. A high error rate indicates large overlap.
#'
#' @param df A data frame containing data for analysis, including metrics and annotations for
#' negative and positive controls.
#' @param metric A string specifying the column name of the metric to use for computing error rate
#' (default: `"log2_normalized_n"`).
#' @param group_cols A character vector specifying the grouping columns (default: `c("depmap_id", "pcr_plate")`).
#' @param negcon A string specifying the `pert_type` value that identifies negative control samples
#' (default: `"ctl_vehicle"`).
#'   These represent baseline or no-treatment conditions and are expected to have lower values.
#' @param poscon A string specifying the `pert_type` value that identifies positive control samples
#' (default: `"trt_poscon"`).
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
#' @param group_cols A character vector specifying the columns to group by, along with `pert_type`
#' (default: `c("depmap_id", "pcr_plate")`).
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
compute_ctl_medians_and_mad <- function(df, group_cols,
                                        negcon = "ctl_vehicle", poscon = "trt_poscon", pseudocount = 20) {
  message("Adding median and MAD values for ", negcon, " and ", poscon, " if it exists.....")
  # Group and compute medians/MADs
  result <- df %>%
    dplyr::filter(pert_type %in% c(negcon, poscon)) %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
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
    )

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
#' @param plate_cell_line_cols A vector of columns describing cell lines.
#' @param sig_cols A vector of columns describing treatment profiles.
#' @return A dataframe.
compute_med_trt_bio_rep = function(df, plate_cell_line_cols, sig_cols) {
  med_trt_bio_reps = df |>
    dplyr::filter(pert_type == "trt_cp") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(unique(c(plate_cell_line_cols, sig_cols))))) |>
    dplyr::summarise(num_trt_bio_reps = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(c(plate_cell_line_cols)))) |>
    dplyr::summarise(med_num_trt_bio_reps = median(num_trt_bio_reps), .groups = "drop")

  return(med_trt_bio_reps)
}

#' Generate cell plate table
#'
#' This function generates a comprehensive QC table for cell_lines +pcr plates by computing and
#' merging various QC metrics, including medians, MADs, error rates, and log fold changes (LFC).
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
generate_plate_cell_table = function(normalized_counts, cell_line_cols, sig_cols,
                                     pseudocount = 20, contains_poscon = TRUE,
                                     pcr_plate_col = "pcr_plate",
                                     pert_plate_col = "pert_plate",
                                     optional_cols = c("project_code", "day"),
                                     poscon = "trt_poscon",
                                     negcon = "ctl_vehicle",
                                     nc_variability_threshold = 1, error_rate_threshold = 0.05,
                                     pc_viability_threshold = 0.25, nc_raw_count_threshold = 40) {
  message("generate_plate_cell_table inputs ...")
  message("  cell_line_cols: ", paste(cell_line_cols, collapse = ", "))
  message("   pcr_plate_col: ", pcr_plate_col)
  message("  pert_plate_col: ", pert_plate_col)
  message("       pert_type: pert_type (hard coded)")
  message("   optional_cols: ", paste(optional_cols, collapse = ", "))
  plate_cell_line_cols = c(cell_line_cols, pcr_plate_col, pert_plate_col)

  # Check that columns used in this module are in normalized counts
  missing_cols = setdiff(c(plate_cell_line_cols, "pert_type"), colnames(normalized_counts))
  if (length(missing_cols) != 0) {
    message("The following columns are missing from normalized_counts: ")
    for (item in missing_cols) {
      message("  ", item)
    }
    stop("Nomalized counts is missing some columns required for the QC tables module.")
  }

  # Pull out optional PCR plate level columns
  optional_meta = normalized_counts |>
    dplyr::distinct(across(all_of(c(pcr_plate_col, pert_plate_col))), across(any_of(optional_cols)))

  # Compute medians and MADs
  message("Calculating medians and MADs ...")
  cell_line_meds_mads = compute_ctl_medians_and_mad(
    df = normalized_counts,
    group_cols = c(plate_cell_line_cols, "pert_type"),
    negcon = negcon,
    poscon = poscon,
    pseudocount = pseudocount
  )

  # Calc median number of bio reps across the treatments for each cell line + plate
  message("Computing median number of bio_reps across the treatments ...")
  med_trt_bio_reps = compute_med_trt_bio_rep(df = normalized_counts,
                                             plate_cell_line_cols = plate_cell_line_cols,
                                             sig_cols = sig_cols)
  if (contains_poscon) {
    # Stats that req poscon
    # Compute error rate - requires pert_type in norm counts
    message("Detected poscons - computing error rates ...")
    error_rates = compute_error_rate(
      df = normalized_counts,
      metric = "log2_normalized_n",
      group_cols = plate_cell_line_cols,
      negcon = negcon,
      poscon = poscon,
      contains_poscon = contains_poscon
    )
    # Join tables together and calculate poscon l2fc
    message("Left joining tables ...")
    plate_cell_table = cell_line_meds_mads |>
      dplyr::left_join(optional_meta, by = c(pcr_plate_col, pert_plate_col)) |>
      dplyr::left_join(med_trt_bio_reps, by = plate_cell_line_cols) |>
      dplyr::left_join(error_rates, by = plate_cell_line_cols) |>
      dplyr::mutate(lfc_trt_poscon = .data[[paste0("median_log_normalized_", poscon)]] -
                      .data[[paste0("median_log_normalized_", negcon)]],
                    viability_trt_poscon = 2^lfc_trt_poscon,
                    qc_pass = error_rate < error_rate_threshold &
                      .data[[paste0("viability_", poscon)]] < pc_viability_threshold &
                      .data[[paste0("median_raw_", negcon)]] > nc_raw_count_threshold &
                      .data[[paste0("mad_log_normalized_", negcon)]] < nc_variability_threshold)

  } else {
    message("No poscon condition detected. Error rate and poscon_l2fc QCs will not be generated.")
    plate_cell_table = cell_line_meds_mads |>
      dplyr::left_join(optional_cols, by = c(pcr_plate_col, pert_plate_col)) |>
      dplyr::left_join(med_trt_bio_reps, by = plate_cell_line_cols) |>
      dplyr::mutate(qc_pass = .data[[paste0("median_raw_", negcon)]] > nc_raw_count_threshold &
                      .data[[paste0("mad_log_normalized_", negcon)]] < nc_variability_threshold)
  }

  # Determine if cell line passes on a pert plate level
  message("Check qc pass rate at pert plate level ...")
  plate_cell_table = plate_cell_table |>
    dplyr::mutate(n_passing_med_num_trt_reps = ifelse(qc_pass, med_num_trt_bio_reps, 0)) |>
    dplyr::group_by(across(all_of(c(cell_line_cols, pert_plate_col)))) |>
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