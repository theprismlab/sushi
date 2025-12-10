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


#' Compute read stats
#'
#' This function calculates summary statistics related to reads and cell line recovery for annotated counts data,
#' combining metadata and performing group-wise computations.
#'
#' @param annotated_counts A data frame containing read data with annotations for barcodes, cell lines, and counts.
#' @param cell_set_meta A data frame containing metadata about cell sets and their expected cell lines.
#' @param group_cols A character vector specifying the grouping columns (default: `c("pcr_plate", "pcr_well")`).
#' @param metric A string indicating the column name of the metric to use for calculations (default: `"n"`).
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line (default: `c("depmap_id", "lua", "pool_id")`).
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered (default: `40`).
#'
#' @return A data frame summarizing read statistics for each group, including total reads, expected reads, control barcode reads,
#' recovered cell lines, and their fractions.
#'
#' @import dplyr
compute_read_stats <- function(annotated_counts, cell_set_meta, unknown_counts, cb_metrics, group_cols = c("pcr_plate", "pcr_well", "pert_type", "pert_plate"),
                               metric = "n", cell_line_cols = c("depmap_id", "pool_id", "lua"), count_threshold = 40,
                               expected_reads_threshold = 0.8, cb_threshold = 40, cb_spearman_threshold = 0.8, cb_mae_threshold = 1) {
  # Compute expected lines from cell_set_meta
  expected_lines <- compute_expected_lines(cell_set_meta, cell_line_cols)

  # Group unknown_counts by group_cols
  unknown_counts <- unknown_counts %>%
    dplyr::left_join(unique(annotated_counts %>% select(pcr_plate, pcr_well, pert_type, pert_plate)),
      by = c("pcr_plate", "pcr_well")
    ) %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      n = sum(.data[[metric]], na.rm = TRUE),
      expected_read = FALSE
    )

  plate_well <- annotated_counts %>%
    dplyr::left_join(expected_lines, by = "cell_set") %>%
    # Append unknown_counts, filling NA for any column not present in unknown_counts
    dplyr::bind_rows(unknown_counts) %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      # Total reads
      n_total_reads = sum(.data[[metric]], na.rm = TRUE),
      # Reads mapping to cell lines
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      # Reads mapping to control barcodes
      n_cb_reads = sum(.data[[metric]][cb_name != ""], na.rm = TRUE),
      # Median reads for control barcodes
      median_cb_reads = median(.data[[metric]][cb_name != ""], na.rm = TRUE),
      # Fraction of reads mapping to cell lines
      fraction_expected_reads = n_expected_reads / n_total_reads,
      # Number of cell lines with coverage above 40 reads
      n_lines_recovered = sum(.data[[metric]] >= count_threshold & (is.na(cb_name) | cb_name == "") & expected_read == TRUE, na.rm = TRUE),
      # Number of expected lines based on metadata
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      # Fraction of cell lines with coverage above count threshold
      fraction_cl_recovered = n_lines_recovered / max(n_expected_lines, na.rm = TRUE),
      # Ratio of control barcode reads to cell line reads
      cb_cl_ratio_well = n_cb_reads / n_expected_reads,
      # Fraction of reads mapping to control barcodes
      fraction_cb_reads = n_cb_reads / n_total_reads
    ) %>%
    dplyr::ungroup()

  plate_pert_type <- plate_well %>%
    dplyr::left_join(cb_metrics %>% select(pcr_plate, pcr_well, cb_mae, cb_spearman), by = c("pcr_plate", "pcr_well")) %>%
    dplyr::filter(median_cb_reads > cb_threshold) %>%
    dplyr::filter(fraction_expected_reads > expected_reads_threshold) %>%
    dplyr::filter(cb_spearman > cb_spearman_threshold & cb_mae < cb_mae_threshold) %>%
    dplyr::group_by(pcr_plate, pert_type) %>%
    dplyr::summarise(
      cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE),
    ) %>%
    dplyr::ungroup()

  # Combine plate_well and plate_pert_type
  result <- plate_well %>%
    dplyr::left_join(plate_pert_type, by = c("pcr_plate", "pert_type"))

  return(result)
}

#' Calculate CB metrics
#'
#' This function calculates control metrics.
#'
#' @param normalized_counts A data frame containing normalized read counts and associated metrics.
#' @param group_cols A character vector specifying the grouping columns
#'                   (default: `c("pcr_plate", "pcr_well", "pert_plate", "cb_intercept")`).
#' @param cb_meta A data frame containing control barcode metadata.
#' @param pseudocount An integer specifying the pseudocount to add to all read counts.
#' @return A data frame containing unique combinations of `group_cols` and the calculated metrics.
#'
#' @import dplyr
#'
calculate_cb_metrics <- function(normalized_counts, cb_meta,
                                 group_cols = c("pcr_plate", "pcr_well", "pert_plate", "cb_intercept"),
                                 pseudocount = 20) {
  # Filter cb_meta by dropping any control barcodes without "well_norm" indicated under "cb_type"
  if ("cb_type" %in% colnames(cb_meta)) {
    dropped_cbs = cb_meta |> dplyr::filter(cb_type != "well_norm")

    if (nrow(dropped_cbs) > 0) {
      print("The following CBs are excluded from normalization.")
      print(dropped_cbs)
      cb_meta = cb_meta |> dplyr::filter(cb_type == "well_norm")
    }
  }

  # Pull out control barcodes used for normalization and calculate some stats
  fit_stats = normalized_counts |>
    dplyr::semi_join(cb_meta, by = c("cb_ladder", "cb_name")) |>
    # To do: Potentially update group_cols to include "log2_pseudovalue"
    dplyr::group_by(across(all_of(group_cols)), across(any_of(c("log2_pseudovalue")))) |>
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

# TABLE GENERATION FUNCTION ----------
#' Generate ID column QC metrics table
#'
#' This function computes and combines various quality control (QC) metrics, such as read statistics, skew,
#' and control barcode metrics, grouped by a specified list of id_cols. For standard sequencing runs, this should always
#' be "pcr_plate" and "pcr_well".
#'
#' @param annotated_counts A data frame containing annotated read counts with metadata.
#' @param normalized_counts A data frame containing normalized read counts and associated metrics/metadata.
#' @param cell_set_meta A data frame containing metadata about cell sets and expected cell lines.
#' @param id_cols_list A character vector specifying the columns to group by for QC metric computations.
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line.
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered (default: `40`).
#'
#' @return A data frame (`id_cols_table`) that combines QC metrics, including read statistics, skew, and control barcode metrics, grouped by the specified id_cols.
#'
#' @import dplyr
generate_id_cols_table <- function(annotated_counts, normalized_counts, unknown_counts, cell_set_meta, cb_meta, id_cols_list, cell_line_cols,
                                   count_threshold = 40, pseudocount = 20, cb_threshold = 40) {
  print(paste0("Computing id_cols QC metrics grouping by ", paste0(id_cols_list, collapse = ","), "....."))

  read_stats_grouping_cols <- c(id_cols_list, "pert_type", "pert_plate")

  cb_metrics <- calculate_cb_metrics(normalized_counts, cb_meta,
                                     group_cols = c(id_cols_list, "pert_plate", "cb_intercept"),
                                     pseudocount = pseudocount)

  read_stats <- compute_read_stats(
    annotated_counts = annotated_counts, unknown_counts = unknown_counts, cb_metrics = cb_metrics, group_cols = read_stats_grouping_cols,
    cell_set_meta = cell_set_meta, metric = "n", cell_line_cols = cell_line_cols,
    count_threshold = count_threshold, cb_threshold = 40
  )

  skew <- compute_skew(annotated_counts, group_cols = c(id_cols_list, "pert_plate"), metric = "n")


  id_cols_table <- read_stats %>%
    dplyr::left_join(skew, by = c(id_cols_list, "pert_plate")) %>%
    dplyr::left_join(cb_metrics, by = c(id_cols_list, "pert_plate")) %>%
    dplyr::left_join(normalized_counts %>% dplyr::select(all_of(c(id_cols_list, "cell_set"))) %>% dplyr::distinct(),
      by = id_cols_list
    ) %>% dplyr::distinct()

  return(id_cols_table)
}

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
    print(paste0("Computing error rate using ", negcon, " and ", poscon, "....."))
    print(paste0("Grouping by ", paste0(group_cols, collapse = ","), "....."))
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
  }
  else {
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
compute_ctl_medians_and_mad <- function(df, group_cols = c("depmap_id", "pcr_plate", "pert_plate"),
                                        negcon = "ctl_vehicle", poscon = "trt_poscon", pseudocount = 20) {
  print(paste0("Adding control median and MAD values for ", negcon, " and ", poscon, " if it exists....."))
  print(paste0("Computing falses sensitivity probability for ", negcon, "....."))
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
    pivot_wider(
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
compute_control_lfc <- function(df, negcon = "ctl_vehicle", poscon = "trt_poscon", grouping_cols = c("depmap_id", "pcr_plate", "pert_plate"), contains_poscon = TRUE) {
  print(paste0("Computing log fold change for ", negcon, " and ", poscon, "....."))
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
  }
  else {
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
    dplyr::group_by(dplyr::pick(tidyselect::all_of(unique(c(cell_line_cols, sig_cols, "pert_plate"))))) |>
    dplyr::summarise(num_trt_bio_reps = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(dplyr::pick(tidyselect::all_of(unique(c(cell_line_cols, "pert_plate"))))) |>
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
generate_cell_plate_table <- function(normalized_counts, filtered_counts, cell_line_cols, sig_cols, pseudocount = 20, contains_poscon = TRUE, poscon = "trt_poscon", negcon = "ctl_vehicle",
                                      nc_variability_threshold = 1, error_rate_threshold = 0.05, pc_viability_threshold = 0.25, nc_raw_count_threshold = 40) {
  cell_line_list <- strsplit(cell_line_cols, ",")[[1]]
  cell_line_plate_grouping <- c(cell_line_list, "pcr_plate", "pert_plate", "project_code", "day") # Define columns to group by
  print(paste0("Computing cell + plate QC metrics grouping by ", paste0(cell_line_plate_grouping, collapse = ","), "....."))

  # Compute control medians and MAD
  medians_and_mad <- compute_ctl_medians_and_mad(
    df = normalized_counts,
    group_cols = cell_line_plate_grouping,
    negcon = negcon,
    poscon = poscon,
    pseudocount = pseudocount
  )

  # Compute error rate
  error_rates <- compute_error_rate(
    df = normalized_counts,
    metric = "log2_normalized_n",
    group_cols = cell_line_plate_grouping,
    negcon = negcon,
    poscon = poscon,
    contains_poscon = contains_poscon
  )

  # Compute poscon LFC
  poscon_lfc <- compute_control_lfc(
    df = medians_and_mad,
    negcon = negcon,
    poscon = poscon,
    grouping_cols = cell_line_plate_grouping,
    contains_poscon = contains_poscon
  )

  # Compute cell line fractions per plate
  cell_line_fractions <- compute_cl_fractions(
    df = filtered_counts,
    grouping_cols = cell_line_plate_grouping
  )

  # Calc median number of bio reps across the treatments for each cell line + plate
  med_trt_bio_reps = compute_med_trt_bio_rep(norm_counts = normalized_counts,
                                             cell_line_cols = cell_line_plate_grouping,
                                             sig_cols = sig_cols)

  # Merge all tables together
  if (contains_poscon) {
    print(paste0("Merging ", paste0(cell_line_plate_grouping, collapse = ","), " QC tables together....."))
    print("medians_and_mad")
    print(colnames(medians_and_mad))
    print("error_rates")
    print(colnames(error_rates))
    print("poscon_lfc")
    print(colnames(poscon_lfc))
    print("cell_line_fractions")
    print(colnames(cell_line_fractions))
    plate_cell_table <- medians_and_mad %>%
      dplyr::left_join(error_rates, by = cell_line_plate_grouping) %>%
      dplyr::left_join(poscon_lfc, by = cell_line_plate_grouping) %>%
      dplyr::left_join(cell_line_fractions, by = cell_line_plate_grouping) %>%
      dplyr::left_join(med_trt_bio_reps, by = cell_line_plate_grouping)
    # QC pass criteria, currently with hardcoded pert_types
    plate_cell_table <- plate_cell_table %>%
      dplyr::mutate(qc_pass = error_rate < error_rate_threshold &
                      viability_trt_poscon < pc_viability_threshold &
                      median_raw_ctl_vehicle > nc_raw_count_threshold &
                      mad_log_normalized_ctl_vehicle < nc_variability_threshold,
                    n_passing_med_num_trt_reps = ifelse(qc_pass, med_num_trt_bio_reps, 0))  %>%
      dplyr::group_by(across(all_of(c(cell_line_list, "pert_plate")))) %>%
      dplyr::mutate(qc_pass_pert_plate = sum(n_passing_med_num_trt_reps) > 1) %>%
      dplyr::ungroup()
    # Add the n_expected_controls values
    plate_cell_table <- plate_cell_table %>%
      dplyr::left_join(
        n_expected_controls,
        by = c("pcr_plate", "pert_plate")
      ) %>%
      dplyr::mutate(
        fraction_expected_poscon = n_replicates_trt_poscon/n_expected_trt_poscon,
        fraction_expected_negcon = n_replicates_ctl_vehicle/n_expected_ctl_vehicle
      )
  }
  else {
    print(paste0("Merging ", paste0(cell_line_plate_grouping, collapse = ","), " QC tables together....."))
    print("medians_and_mad")
    print(colnames(medians_and_mad))
    print("cell_line_fractions")
    print(colnames(cell_line_fractions))
    plate_cell_table <- medians_and_mad %>%
      dplyr::left_join(cell_line_fractions, by = cell_line_plate_grouping) %>%
      dplyr::left_join(med_trt_bio_reps, by = cell_line_plate_grouping)
    # QC pass criteria, currently with hardcoded pert_types
    plate_cell_table <- plate_cell_table %>%
      dplyr::mutate(qc_pass = median_raw_ctl_vehicle > nc_raw_count_threshold &
                      mad_log_normalized_ctl_vehicle < nc_variability_threshold,
                    n_passing_med_num_trt_reps = ifelse(qc_pass, med_num_trt_bio_reps, 0))  %>%
      dplyr::group_by(across(all_of(c(cell_line_list, "pert_plate")))) %>%
      dplyr::mutate(qc_pass_pert_plate = sum(n_passing_med_num_trt_reps) > 1) %>%
      dplyr::ungroup()
    # Add the n_expected_controls values
    plate_cell_table <- plate_cell_table %>%
      dplyr::left_join(
        n_expected_controls,
        by = c("pcr_plate", "pert_plate")
      ) %>%
      dplyr::mutate(
        fraction_expected_negcon = n_replicates_ctl_vehicle/n_expected_ctl_vehicle
      )
  }

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
  qc_table <- plate_cell_table %>%
    mutate(qc_flag = case_when(
      mad_log_normalized_ctl_vehicle > nc_variability_threshold ~ "nc_variability",
      error_rate > error_rate_threshold ~ "error_rate",
      viability_trt_poscon > pc_viability_threshold ~ "pc_viability",
      median_raw_ctl_vehicle < log(nc_raw_count_threshold) ~ "nc_raw_count",
      TRUE ~ NA_character_
    ))
  }
  else {
    qc_table <- plate_cell_table %>%
      mutate(qc_flag = case_when(
        mad_log_normalized_ctl_vehicle > nc_variability_threshold ~ "nc_variability",
        median_raw_ctl_vehicle < log(nc_raw_count_threshold) ~ "nc_raw_count",
        TRUE ~ NA_character_
      ))
  }
  return(qc_table)
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
                             contamination_threshold = contamination_threshold,
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
generate_pool_well_qc_table <- function(normalized_counts, pool_well_delta_threshold = 5, pool_well_fraction_threshold = 0.4) {
  ### POOL_WELL_OUTLIERS
  ## Flag pool/well combinations based on the fraction of cell lines in a pool + well that are some distance from the pool + well median.
  ## This is in control wells only
  pool_well_outliers <- normalized_counts %>%
    # Consider only cell line barcodes
    filter(is.na(cb_name)) %>%
    # Get the median value of the pool in each well
    group_by(pcr_plate, pcr_well, pert_type, pool_id) %>%
    mutate(
      pool_well_median = median(log2_normalized_n, na.rm = TRUE),
      n_cell_lines = n_distinct(paste(lua, depmap_id, cell_set, sep = "_"))
    ) %>%
    mutate(
      delta_from_pool_well_median = abs(log2_normalized_n - pool_well_median)
    ) %>%
    group_by(pcr_plate, pcr_well, pert_type, pool_id) %>%
    reframe(
      n_outliers = sum(delta_from_pool_well_median > pool_well_delta_threshold, na.rm = TRUE),
      n_cell_lines = max(n_cell_lines)
    ) %>%
    ungroup() %>%
    mutate(fraction_outliers = n_outliers / n_cell_lines,
      qc_flag = if_else(fraction_outliers > pool_well_fraction_threshold, "pool_well_outliers", NA_character_))
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
  }
  else {
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

# Variance decomposition table
compute_variance_decomposition <- function(normalized_counts, metric = 'n', negcon = "ctl_vehicle",
                                           cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
                                           id_cols = c("pcr_plate","pcr_well")) {
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

# CONTAMINATION QC TABLES

compute_contamination_qc_tables <- function(prism_barcode_counts,
                                            unknown_barcode_counts,
                                            cell_set_and_pool_meta,
                                            cell_line_meta,
                                            cb_meta,
                                            sample_meta) {
  # --- 0. Ensure all required columns are present ---
  prism_barcode_counts_all <- prism_barcode_counts %>%
    dplyr::left_join(sample_meta %>% select(pcr_plate, pcr_well, pert_plate),
         by = c("pcr_plate", "pcr_well"))

  unknown_barcode_counts_all <- unknown_barcode_counts %>%
      dplyr::left_join(sample_meta %>% select(pcr_plate, pcr_well, pert_plate),
           by = c("pcr_plate", "pcr_well"))

  # --- 1. Create total counts df with known and unknown barcodes ---
  total_counts <- bind_rows(prism_barcode_counts_all, unknown_barcode_counts_all)

  # --- 2. Compute total counts for each well ---
  total_counts_by_well <- total_counts %>%
    group_by(pcr_plate, pcr_well, pert_plate) %>%
    summarise(well_count = sum(n), .groups = "drop")

  # --- 3. Get a list of expected reads ---
  expected_cell_lines <- cell_set_and_pool_meta %>%
    left_join(cell_line_meta, by = "depmap_id") %>%
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
  total_counts_by_read_type <- total_counts_with_read_type_annotated %>% select(c("pcr_plate","pert_plate","pcr_well","forward_read_barcode","read_type","n"))

  return(list(fraction_read_type_by_well = fraction_read_type_by_well,
              fraction_read_type_by_plate = fraction_read_type_by_plate,
              total_counts_by_read_type = total_counts_by_read_type))
}

