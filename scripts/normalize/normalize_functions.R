#' get_valid_norm_cbs
#'
#' Extracts the control barcodes from the input filtered counts table and
#' identifies control barcodes that can be used for normalization.
#'
#' This function flags control barcodes that are not present in the CB meta table,
#' undetected in sequencing, and have high MAD in the negative control wells of the PCR plate.
#' It also flags PCR wells that do not have enough usable control barcodes (<= 4).
#'
#' @param filtered_counts Data table of filtered counts.
#' @param CB_meta Data table of the control barcode metadata.
#' @param id_cols Vector of columns that uniquely identify each PCR well.
#' @param negcon_type String identifying the negative controls in the pert_type column.
#' @param cb_mad_cutoff Numeric maximum MAD value for the control barcodes.
#' @param req_negcon_reps Integer number of negative control replicates required for the MAD filter.
#' @return Data frame of control barcodes with column indicating status or failure mode.
#' @import tidyverse
flag_control_bcs = function(filtered_counts, CB_meta, id_cols, negcon_type,
                            cb_mad_cutoff = 1, req_negcon_reps = 6) {
  # Drop any control barcodes in CB_meta NOT marked with "well_norm".
  if ("cb_type" %in% colnames(CB_meta)) {
    dropped_cbs = CB_meta |> dplyr::filter(cb_type != "well_norm")

    if (nrow(dropped_cbs) > 0) {
      message("The following CBs in CB_meta are excluded from normalization.")
      print(dropped_cbs)
      CB_meta = CB_meta |> dplyr::filter(cb_type == "well_norm")

      if (nrow(CB_meta) == 0) {
        message("There are no CBs in CB_meta that can be used for normalization.")
        stop("CB_meta needs a control barcode ladder marked with 'well_norm' in the 'cb_type' column.")
      }
    }
  }

  # Create a CB flag column
  # This will later be merged onto a df of all CBs.
  small_CB_meta = CB_meta |>
    dplyr::distinct(cb_ladder, cb_name) |>
    dplyr::mutate(keep_cb = "Yes")

  # Flag control barcodes that are:
  # 1. NOT present in CB_meta for the screen
  # 2. undetected in sequencing
  cbs_annots = filtered_counts |>
    dplyr::filter(!is.na(cb_name)) |>
    dplyr::distinct(dplyr::across(tidyselect::all_of(c(id_cols, "pert_type", "cb_ladder", "cb_name", "n")))) |>
    dplyr::left_join(small_CB_meta, by = c("cb_ladder", "cb_name")) |>
    dplyr::mutate(keep_cb = dplyr::case_when(
      is.na(keep_cb) ~ "Not in CB meta", # Flag CBs not in CB_meta
      n == 0 ~ "Undetected CB", # Flag undetected CBs
      .default = "Yes")
    )
  # Report if any CBs were flagged
  total_cbs = nrow(cbs_annots)
  not_in_meta_cbs = nrow(cbs_annots[keep_cb == "Not in CB meta"])
  undetected_cbs = nrow(cbs_annots[keep_cb == "Undetected CB"])
  message(sprintf("%d (%.2f%%) control barcodes were not found in the CB meta.",
                  not_in_meta_cbs, not_in_meta_cbs / total_cbs))
  message(sprintf("%d (%.2f%%) control barcodes were not detected in sequencing.",
                  undetected_cbs, undetected_cbs / total_cbs))

  # Flag control barcodes that:
  # 3. have high MAD (variability) in the negative controls of a PCR plate.
  # If there are negative controls, then calculate MADs of the CBs across the negative controls
  # and flag CBs with high MADs on PCR plates with enough negative controls.
  if (negcon_type %in% unique(cbs_annots$pert_type)) {
    # Calculate negcon MADs of CBs for normalization
    cb_mad = get_cb_mad(cbs_annots[pert_type == negcon_type & keep_cb != "Not in CB meta"],
                        id_cols = id_cols,
                        cb_mad_cutoff = cb_mad_cutoff)
    high_mad_cbs = cb_mad |>
      dplyr::filter(keep_cb == FALSE, num_reps >= req_negcon_reps) |>
      dplyr::rename(cb_mad_flag = keep_cb)
    # TO DO: consider different control conditions - different days or different conditions.

    # Flag barcodes with high MAD if there are any
    if (nrow(high_mad_cbs) > 0) {
      cbs_annots = cbs_annots |>
        dplyr::left_join(high_mad_cbs, by = c("pcr_plate", "cb_name"), suffix = c("", ".y")) |>
        dplyr::select(!tidyselect::ends_with(".y")) |>
        dplyr::mutate(keep_cb = dplyr::case_when(cb_mad_flag == FALSE & keep_cb == "Yes" ~ "High negcon MAD",
                                                 .default = keep_cb))
    }
    # Report the number of CBs flagged and the number of CBs that were not filtered due to low replicate count
    message(sprintf("%d (%.2f%%) control barcodes have MADs greater than %.2f.",
                    nrow(high_mad_cbs), nrow(high_mad_cbs) / total_cbs, cb_mad_cutoff))
    message(sprintf("The MAD filter did not apply to %d CBs where the number of replicates is below %d.",
                    nrow(cb_mad |> dplyr::filter(keep_cb == FALSE, num_reps < req_negcon_reps)),
                    req_negcon_reps))
  } else {
    message("No negative controls detected. Skipping the CB MAD filter.")
  }

  # Flag control barcodes that:
  # 4. Flag PCR wells without enough unflagged control barcodes for normalization.
  # In a PCR well if there are 4 or fewer CBs tagged with "Yes" in keep_cb,
  # those "Yes" notes are converted into "Not enough valid CBs".
  cbs_annots = cbs_annots |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(id_cols))) |>
    dplyr::mutate(keep_cb = ifelse(sum(keep_cb == "Yes") <= 4 & keep_cb == "Yes",
                                   "Not enough valid CBs", keep_cb)) |>
    dplyr::ungroup()

  # Print out the number of PRC wells that can be normalized
  num_pcr_wells = filtered_counts |> dplyr::distinct(dplyr::across(tidyselect::all_of(id_cols))) |> nrow()
  num_passing_wells = cbs_annots |>
    dplyr::filter(keep_cb == "Yes") |>
    dplyr::distinct(dplyr::across(tidyselect::all_of(id_cols))) |>
    nrow()
  message(sprintf("Out of %d PCR wells, %d wells contain enough CBs for normalization.",
                  num_pcr_wells, num_passing_wells))

  # Throw an error if there are no valid PCR wells to normalize
  if (num_passing_wells == 0) {
    flag_summary = cbs_annots |>
      dplyr::group_by(keep_cb) |>
      dplyr::summarise(num_cbs = dplyr::n(), .groups = "drop")
    print(flag_summary)
    stop("No valid PCR wells identified for normalization.")
  }

  return(cbs_annots)
}

#' normalize
#'
#' Normalize read counts using valid control barcodes.
#'
#' @param X Data table of filtered read counts.
#' @param valid_cbs Data frame of the control barcodes with a column "keep_cb".
#' @param id_cols Vector of columns that uniquely identify each PCR well.
#' @param pseudocount Integer to be added to all counts so that logs can be taken.
#' @return Data frame of normalized counts.
#' @import tidyverse
normalize = function(X, cb_annots, id_cols, pseudocount = 0) {
  # Validation: Check that id_cols are present in the dataframe
  if (!validate_columns_exist(id_cols, X)) {
    stop("One or more id_cols (printed above) is NOT present in filtered counts.")
  }

  # Drop cell lines that appear in more than one pool.
  # These lines were identified in the filtered counts module and have a specific string
  # in the pool_id column.
  if ("pool_id" %in% colnames(X)) {
    X = X |> dplyr::filter(!grepl("_+_", pool_id, fixed = TRUE))
  }

  # Add pseudocount and join control barcode annotations
  X = X |>
    dplyr::mutate(log2_n = log2(n + pseudocount)) |>
    dplyr::left_join(cb_annots, by = c(id_cols, "cb_name"), suffix = c("", ".y")) |>
    dplyr::select(!tidyselect::ends_with(".y")) # Drop any duplicate columns from cb_annots

  # Calculate fit intercept for valid profiles using median intercept
  fit_intercepts = X |>
    dplyr::filter(keep_cb == "Yes") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(c(id_cols, "cb_log2_dose")))) |>
    dplyr::summarize(dose_intercept = mean(cb_log2_dose) - mean(log2_n), .groups = "drop") |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(id_cols))) |>
    dplyr::summarize(cb_intercept = median(dose_intercept), .groups = "drop")

  # Normalize filtered read counts and note flagged CBs
  normalized = dplyr::inner_join(X, fit_intercepts, by = id_cols) |>
    dplyr::mutate(log2_normalized_n = log2_n + cb_intercept,
                  cb_ladder = dplyr::case_when(is.na(keep_cb) ~ cb_ladder,
                                               keep_cb == "Yes" ~ cb_ladder,
                                               keep_cb != "Yes" ~ paste(cb_ladder, keep_cb, sep = " - "))) |>
    dplyr::select(-log2_n, -keep_cb)

  return(normalized)
}

#' get_cb_mad
#'
#' Calculates MAD for each control barcodes across all negative control wells of a plate
#'
#' @param cb_reads Dataframe of control barcode read counts.
#' @param id_cols Vector of columns that uniquely identify each PCR well.
#' @param cb_mad_cutoff Maximum MAD value for a control barcode.
get_cb_mad = function(cb_reads, id_cols, cb_mad_cutoff = 1) {
  cb_mad = cb_reads |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(id_cols))) |>
    dplyr::mutate(cb_frac = (n + 1) / sum(n)) |>
    dplyr::group_by(pcr_plate, cb_ladder, cb_name) |>
    dplyr::summarise(mad_log2_cb_frac = mad(log2(cb_frac)),
                     num_reps = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(keep_cb = ifelse(mad_log2_cb_frac < cb_mad_cutoff, TRUE, FALSE))

  return(cb_mad)
}

#' add_pseudovalues
#'
#' After normalization without a pseudocount, compute a pseudovalue using the negative controls.
#' Add the pseudovalue to normalized counts.
#'
#' @param norm_counts Dataframe of normalized counts.
#' @param negcon_cols Vector of columns names in norm_counts that describe a negative control condition.
#' @param read_detection_limit Read count detection limit.
#' @param negcon_type String in column pert_type of norm_counts that indicates negative control samples.
#' @returns Dataframe with normalized count adjusted with a pseudovalue.
#' @import tidyverse
add_pseudovalue = function(norm_counts, negcon_cols, read_detection_limit = 10, negcon_type = "ctl_vehicle") {
  # Pseudovalue is calculated over negative control groups.
  # Check that all of those columns are present.
  missing_cols = setdiff(negcon_cols, colnames(norm_counts))
  if (length(missing_cols) > 0) {
    stop("normalize - add_pseudovalue: The following column(s) are missng from normalized_counts: ",
         paste(missing_cols, collapse = ","))
  }

  # Calculate the pseudovalue over each control group.
  pv_per_group = norm_counts |> dplyr::filter(pert_type == negcon_type) |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(negcon_cols))) |>
    dplyr::summarise(log2_pseudovalue = median(log2(read_detection_limit) + cb_intercept), .groups = "drop")

  # Join pseudovalues to the main df and add them to all rows/samples.
  norm_counts = norm_counts |>
    dplyr::left_join(pv_per_group, by = negcon_cols) |>
    dplyr::mutate(log2_normalized_n = log2(2^log2_normalized_n + 2^log2_pseudovalue))

  return(norm_counts)
}