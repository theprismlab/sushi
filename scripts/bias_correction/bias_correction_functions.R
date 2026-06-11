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



# Load required libraries
library(dplyr)
library(rlang)

#' Variance-Adaptive Shrinkage Correction
#' 
#' @param df Dataframe in long format
#' @param grouping_vars Character vector of columns defining experimental units (e.g., c("drug_id"))
#' @param corrected_val String name of the column with OLS corrected values
#' @param val String name of the column with raw uncorrected values
#' @param batch_vars Character vector of columns defining discrete artifacts (e.g., c("a", "b"))
#' @param shrunk_col_name String name for the output column (defaults to "shrunk_val")
#' @return The original dataframe with `lambda` and `shrunk_val` appended.
apply_shrinkage_correction_heavy <- function(df, 
                                       grouping_vars = c("bio_rep", "CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "day"),
                                       corrected_val = "l2fc" , 
                                       val = "l2fc_uncorrected", 
                                       batch_vars = c("cell_set", "growth_condition"),
                                       shrunk_col_name = "lfc_adaptive") {
  
  # 1. Convert strings to symbols for tidy evaluation
  grouping_syms <- syms(grouping_vars)
  batch_syms <- syms(batch_vars)
  corrected_sym <- sym(corrected_val)
  val_sym <- sym(val)
  
  # 2. Define the ultra-fast, fail-safe lambda calculator
  
  # calc_lambda <- function(v, g) {
  #   # Filter NAs to prevent math errors
  #   valid <- !is.na(v) & !is.na(g)
  #   v <- v[valid]
  #   g <- g[valid]
  #   
  #   # Safety: need at least 2 groups and 3 points to calculate variance
  #   if (length(unique(g)) < 2 || length(v) < 3) return(1)
  #   
  #   # Step A: Calculate absolute deviations from group medians (Levene's Z)
  #   group_medians <- ave(v, g, FUN = median)
  #   Z <- abs(v - group_medians)
  #   
  #   # Step B: Calculate Sums of Squares directly
  #   grand_mean_Z <- mean(Z)
  #   group_means_Z <- ave(Z, g, FUN = mean)
  #   
  #   SS_total <- sum((Z - grand_mean_Z)^2)
  #   SS_between <- sum((group_means_Z - grand_mean_Z)^2)
  #   
  #   # Safety: If there is literally zero variance in the data
  #   if (SS_total == 0) return(1)
  #   
  #   # Step C: Compute Eta-squared and map to Lambda
  #   eta_squared <- SS_between / SS_total
  #   lambda <- 1 - sqrt(eta_squared)
  #   
  #   # Clamp between 0 and 1 to prevent floating point weirdness
  #   return(max(0, min(1, lambda)))
  # }

  calc_lambda <- function(v, g) {
    valid <- !is.na(v) & !is.na(g)
    v <- v[valid]
    g <- g[valid]
    
    k <- length(unique(g))
    N <- length(v)
    
    if (k < 2 || N < 3) return(1)
    
    group_medians <- ave(v, g, FUN = median)
    Z <- abs(v - group_medians)
    
    grand_mean_Z <- mean(Z)
    group_means_Z <- ave(Z, g, FUN = mean)
    
    SS_total <- sum((Z - grand_mean_Z)^2)
    SS_between <- sum((group_means_Z - grand_mean_Z)^2)
    SS_within <- SS_total - SS_between
    
    if (SS_within < 1e-8) return(0) # Perfect variance separation (extreme ADC)
    
    # Calculate the Levene's W statistic
    W_stat <- (SS_between / (k - 1)) / (SS_within / (N - k))
    
    # Anchor: Find the 99.9th percentile for these specific degrees of freedom
    W_crit <- qf(0.999, df1 = k - 1, df2 = N - k)
    
    # Define the penalty constant so that at W_crit, the penalty is significant
    c_penalty <- 1 / W_crit
    
    # Only penalize W statistics greater than 1 (the expected null value)
    W_adjusted <- max(0, W_stat - 1)
    
    # Calculate lambda (decay function)
    lambda <- 1 / (1 + c_penalty * W_adjusted)
    
    return(max(0, min(1, lambda)))
  }
  
  
  
   
  # 3. Execute the pipeline
  result <- df %>%
    # Create a single grouping variable for the interaction of a and b
    dplyr::mutate(..batch_interaction.. = paste(!!!batch_syms, sep = "_")) %>%
    # Group by drug/compound
    dplyr::group_by(!!!grouping_syms) %>%
    dplyr::mutate(
      lambda = calc_lambda(!!val_sym, ..batch_interaction..),
      # Apply the convex combination
      !!sym(shrunk_col_name) := (!!val_sym) * (1 - lambda) + (!!corrected_sym) * lambda
    ) %>%
    dplyr::ungroup() %>%
    # Clean up the temporary interaction column
    dplyr::select(-..batch_interaction..)
  
  return(result)
}





#' Variance-Adaptive Shrinkage Correction
#' 
#' @param df Dataframe in long format
#' @param grouping_vars Character vector of columns defining experimental units (e.g., c("drug_id"))
#' @param corrected_val String name of the column with OLS corrected values
#' @param val String name of the column with raw uncorrected values
#' @param batch_vars Character vector of columns defining discrete artifacts (e.g., c("a", "b"))
#' @param shrunk_col_name String name for the output column (defaults to "shrunk_val")
#' @return The original dataframe with `lambda` and `shrunk_val` appended.
apply_shrinkage_correction_mild <- function(df, 
                                       grouping_vars = c("bio_rep", "CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "day"),
                                       corrected_val = "l2fc" , 
                                       val = "l2fc_uncorrected", 
                                       batch_vars = c("cell_set", "growth_condition"),
                                       shrunk_col_name = "lfc_adaptive") {
  
  # 1. Convert strings to symbols for tidy evaluation
  grouping_syms <- syms(grouping_vars)
  batch_syms <- syms(batch_vars)
  corrected_sym <- sym(corrected_val)
  val_sym <- sym(val)
  
  # 2. Define the ultra-fast, fail-safe lambda calculator

  calc_lambda <- function(v, g) {
    # Filter NAs to prevent math errors
    valid <- !is.na(v) & !is.na(g)
    v <- v[valid]
    g <- g[valid]

    # Safety: need at least 2 groups and 3 points to calculate variance
    if (length(unique(g)) < 2 || length(v) < 3) return(1)

    # Step A: Calculate absolute deviations from group medians (Levene's Z)
    group_medians <- ave(v, g, FUN = median)
    Z <- abs(v - group_medians)

    # Step B: Calculate Sums of Squares directly
    grand_mean_Z <- mean(Z)
    group_means_Z <- ave(Z, g, FUN = mean)

    SS_total <- sum((Z - grand_mean_Z)^2)
    SS_between <- sum((group_means_Z - grand_mean_Z)^2)

    # Safety: If there is literally zero variance in the data
    if (SS_total == 0) return(1)

    # Step C: Compute Eta-squared and map to Lambda
    eta_squared <- SS_between / SS_total
    lambda <- 1 - sqrt(eta_squared)

    # Clamp between 0 and 1 to prevent floating point weirdness
    return(max(0, min(1, lambda)))
  }
  

  
  
  # 3. Execute the pipeline
  result <- df %>%
    # Create a single grouping variable for the interaction of a and b
    dplyr::mutate(..batch_interaction.. = paste(!!!batch_syms, sep = "_")) %>%
    # Group by drug/compound
    dplyr::group_by(!!!grouping_syms) %>%
    dplyr::mutate(
      lambda = calc_lambda(!!val_sym, ..batch_interaction..),
      # Apply the convex combination
      !!sym(shrunk_col_name) := (!!val_sym) * (1 - lambda) + (!!corrected_sym) * lambda
    ) %>%
    dplyr::ungroup() %>%
    # Clean up the temporary interaction column
    dplyr::select(-..batch_interaction..)
  
  return(result)
}
