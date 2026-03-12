# Collection of functions used by the synergy module

#' restructure_l2fc
#'
#' In the context of a CPS, this function restructures a l2fc file by creating additional
#' columns containing the l2fc values for the single agents of a combination.
#'
#' This function splits the l2fc data.table into rows with single agents and rows with combinations using
#' `combination_col`. The single agent rows are joined onto the combination rows -
#' once for each member of the combination. The new l2fc columns created by the joins
#' are named following `new_l2fc_col_names`.
#'
#' @import data.table
#' @param cps_l2fc
#' @param pert_cols_list A list of vectors where each item is a vector of column names describing a perturbation.
#' @param join_cols A vector of additional column names used for joining.
#' @param l2fc_col A string name of the column containing the l2fc values.
#' @param new_l2fc_col_names A vector of l2fc column names to be created by this function.
#' @param combination_col Name of the column that identifies a treatment as a combination.
#' @return A data.table
restructure_l2fc = function(cps_l2fc, pert_cols_list, join_cols,
                            l2fc_col = "median_l2fc",
                            new_l2fc_col_names = c("pert_l2fc", "pert2_l2fc", "combo_l2fc"),
                            combination_col = "is_combination") {

  # Check that combination_col is a column in cps_l2fc
  if (!combination_col %in% colnames(cps_l2fc)) {
    stop(paste0("The combination_col ", combination_col, " is not in the provided data.table."))
  }

  # Check that vectors in pert_cols_list are of the same size
  if (length(unique(lengths(pert_cols_list))) !=  1) {
    print(pert_cols_list)
    stop("Number of columns describing each perturbation does not match.")
  }

  # Slice dataframe into single agents and combos
  singles_df = cps_l2fc[get(combination_col) == FALSE, ]
  # In combos_df rename l2fc_col to last item of new_l2fc_col_names
  num_names = length(new_l2fc_col_names)
  combos_df = cps_l2fc[get(combination_col) == TRUE,
                       data.table::setnames(.SD, l2fc_col, new_l2fc_col_names[num_names])]

  # For each single agent, add its l2fc column onto combos_df
  # Loop over new_col_names ignoring last item
  for (idx in seq_along(new_l2fc_col_names[-num_names])) {
    # The first iterations will not be a cross join, but subsequent joins will be
    cross_join_cols = setNames(pert_cols_list[[1]], pert_cols_list[[idx]])
    combos_df[singles_df, new_l2fc_col_names[idx] := get(l2fc_col), on = c(join_cols, cross_join_cols)]
  }

  # Reorder columns - move combo l2fc to after single l2fcs
  data.table::setcolorder(combos_df, neworder = new_l2fc_col_names[num_names],
                          after = new_l2fc_col_names[num_names - 1])

  return(combos_df)
}

#' Calculate synergy scores
#'
#' In the context of a CPS, this function calculates synergy scores from l2fc values.
#'
#' From the restructure l2fc data.table, l2fc columns are converted to viabilities which are then used
#' to calculate HSA, Bliss, and synergy scores. Only the synergy related columns are reported.
#' This function relies on the `data.table` package for speed and thus will modify the `restructured_l2fc` input.
#'
#' @import data.table
#' @param restructure_l2fc
#' @param l2fc_cols Vector of column names containing l2fc values.
#' @param viab_cap Upper bound for the viability. Defaults to 1.
#' @return A data.table
calculate_synergy = function(restructured_l2fc, l2fc_cols, viab_cap = 1) {

  # Calculate synergy score from HSA and Bliss
  restructured_l2fc[, hsa := pmin(2^pmin(get(l2fc_cols[1]), get(l2fc_cols[2])), viab_cap)]
  restructured_l2fc[, bliss := pmin(2^get(l2fc_cols[1]), viab_cap) * pmin(2^get(l2fc_cols[2]), viab_cap)]
  restructured_l2fc[, combo_viab := pmin(2^get(l2fc_cols[3]), viab_cap)]
  restructured_l2fc[, synergy := data.table::fcase(combo_viab < bliss, bliss - combo_viab,
                                                   combo_viab >= bliss & combo_viab <= hsa, 0,
                                                   combo_viab > hsa, hsa - combo_viab,
                                                   default = NA)]

  restructured_l2fc[, combo_viab := NULL] # drop combo_viab column

  return(restructured_l2fc)
}

create_dmso_synergy_hdf5 = function(dmso_l2fc, group_name_cols, path,
                                    n_samples, size = 3, replace = FALSE, seed = 2) {
  # Convert dmso_l2fc to data.table if it is not
  if (data.table::is.data.table(dmso_l2fc) == FALSE) {
    message("Converting dmso_l2fc to data.table object.")
    dmso_l2fc = data.table::as.data.table(dmso_l2fc)
  }

  # Remove exisitng hdf5 file of the same name
  if (file.exists(path)) {
    message("Removing existing h5 file at ", path)
    file.remove(path)
  }

  message("Creating new h5 file at ", path)
  rhdf5::h5createFile(path)

  # Create unique group names for each plate + cell line for hdf5 hierarchy
  dmso_l2fc[, group_name := do.call(paste, c(.SD, sep = "__")), .SDcols = group_name_cols]

  # Loop through each unique name
  unique_group_names = unique(dmso_l2fc$group_name)
  for (i in unique_group_names) {
    subset = dmso_l2fc[group_name == i, ]

    # Check resampling maximum
    num_pick_combinations = base::choose(nrow(subset), size)
    if (num_pick_combinations < n_samples) {
      message(nrow(subset), " choose ", size, " = ", num_pick_combinations)
      message("WARNING: Cannot sample up to ", n_samples, ". Using ", num_pick_combinations, " instead.")
      resample_n = num_pick_combinations
    } else {
      resample_n = n_samples
    }

    # Resample to n_samples and create pert1, pert2, combo
    resampled_l2fc = median_resample(x = subset$l2fc, n_samples = resample_n, size = size,
                                     replace = replace, seed = seed)
    mock_values = data.table(mock1 = resampled_l2fc)
    set.seed(seed + 1L)
    mock_values[, mock2 := sample(mock1, size = resample_n, replace = TRUE)]
    set.seed(seed + 2L)
    mock_values[, mock3 := sample(mock1, size = resample_n, replace = TRUE)]

    # Calculate synergy - function occurs in place
    calculate_synergy(restructured_l2fc = mock_values, l2fc_cols = colnames(mock_values))

    # Write to hdf5
    rhdf5::h5write(obj = mock_values$synergy, file = path, name = i)
  }

  print("Closing h5 file.")
  rhdf5::h5closeAll()
}

# Sample DMSO l2fc values
#' Median resample
#'
#' Resample up a larger sample size by picking a certain number of samples and taking the median.
#'
#' This uses `base::sample` and `base::replicate` to resample.
#'
#' @param x Vector of the initial sample values.
#' @param n_samples New sample size to resample to.
#' @param size Number of initial samples to pick to form a group.
#' @param seed Random seed
#' @param replace Sampling with or without replacement for `size`.
#' @param prob Optional vector of probability weights for `base::sample`.
median_resample = function(x, n_samples, size = 3, seed = 2, replace = FALSE, prob = NULL) {
  # Reset seed and resample
  base::set.seed(seed)
  resampled_values = base::replicate(n_samples, median(base::sample(x, size, replace = replace, prob = prob)))

  return(resampled_values)
}

#' Calculate pvalue for a synergy score
#'
#' Calculates a pvalue for a given synergy score using a vector of values representing a null distribution.
#'
#' @import rhdf5
#' @param group_name String name in the `hf_file` hierarchy.
#' @param synergy_value Synergy value
#' @param h5_file hdf5 file containing values for null distributions.
get_pvalue = function(group_name, synergy_value, h5_file) {
  ecdf_obj = stats::ecdf(rhdf5::H5Dread(rhdf5::H5Dopen(h5_file, name = group_name)))
  cdf_value = smooth_pvalue(ecdf_obj(synergy_value), environment(ecdf_obj)$nobs)
  pvalue = 2 * pmin(1 - cdf_value, cdf_value)

  return(pvalue)
}

#' Smooth pvalue
#'
#' Use Laplace smoothing on a pvalue to prevent zero values.
#'
#' @param pval_naive Initial pvalue.
#' @param n_samples Size of the null distribution.
smooth_pvalue = function(pval_naive, n_samples) {
  return((pval_naive * n_samples + 1) / (n_samples + 2))
}
