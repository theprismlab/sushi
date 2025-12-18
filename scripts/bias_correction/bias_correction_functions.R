#' Apply l2fc correction
#'
#' Corrects l2fc change values by reducing variance from growth patterns,
#' cell sets, and abundance in the negative controls. This assumes that cell sets
#' are on different PCR plates, so regressing this out may also be the same as regressing
#' out PCR plate.
#'
#' @param df dataframe
#' @param raw_l2fc_col String name of a column in df containing the l2fc values.
#' @param growth_pattern_col String name of a column in df containing the growth annotations.
#' @param negcon_norm_col String name of a column in df containing negcon_log2_norm_n.
#' @param cell_set_col String name of a column in df containing cell set information.
apply_bias_correction = function(df, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_condition",
                                 negcon_norm_col = "negcon_log2_norm_n",
                                 cell_set_col = "cell_set") {

  # Filter out infinite/NA rows and set factors
  centered_df = df |> dplyr::filter(is.finite(.data[[raw_l2fc_col]]), is.finite(.data[[negcon_norm_col]]),
                                    !is.na(.data[[growth_pattern_col]]), !is.na(.data[[cell_set_col]]))
  centered_df[[growth_pattern_col]] = as.factor(centered_df[[growth_pattern_col]])
  centered_df[[cell_set_col]] = as.factor(centered_df[[cell_set_col]])
  centered_df$centered_l2fc = centered_df[[raw_l2fc_col]] - mean(centered_df[[raw_l2fc_col]])

  # Identify number of growth patterns
  num_growth_patterns = length(unique(centered_df[[growth_pattern_col]]))
  num_cell_sets = length(unique(centered_df[[cell_set_col]]))

  model_formula = create_model_formula(y = "centered_l2fc",
                                       a = growth_pattern_col, count_a = num_growth_patterns,
                                       b = negcon_norm_col,
                                       c = cell_set_col, count_c = num_cell_sets)
  message("Correction model: ", model_formula)

  # Fit model and adjust l2fcs
  fit = lm(formula = as.formula(model_formula), data = centered_df)
  centered_df$l2fc_uncorrected = centered_df$l2fc
  centered_df$l2fc = fit$residuals + mean(centered_df$l2fc)

  # Drop columns created by this function and return output
  return(centered_df |> dplyr::select(-centered_l2fc))
}

#' Create l2fc correction model
#'
#' Model is of the form y = a + b * c, but will change if there is just one group of a or b
create_model_formula = function(y, a, b, c, count_a, count_c) {
  if (count_a > 1) {
    if (count_c > 1) {
      model_formula = sprintf("%s ~ %s + %s * %s", y, a, b, c)
    } else {
      model_formula = sprintf("%s ~ %s + %s", y, a, b)
    }
  } else {
    if (count_c > 1) {
      model_formula = sprintf("%s ~ %s * %s", y, b, c)
    } else {
      model_formula = sprintf("%s ~ %s", y, b)
    }
  }
  return(model_formula)
}