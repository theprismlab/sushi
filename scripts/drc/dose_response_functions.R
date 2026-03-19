#' Computing area under a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param inflection : inflection point (EC50)
#' @param slope : Hill slope ( > 0 for decreasing curves)
#' @param minimum_dose : Minimum dose
#' @param maximum_dose : Maximum dose
#'
#' @return auc value of the dose-response function (x: log2(dose), y: response) scaled to the dose range.
#' @export
#'
#' @examples
compute_auc <- function(LL, UL, inflection, slope, minimum_dose, maximum_dose) {
  f1 = function(x) pmax(pmin((UL + (LL - UL)/(1 + (2^x/inflection)^slope)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(minimum_dose), log2(maximum_dose))$value/(log2(maximum_dose/minimum_dose)),
                  error = function(e) {print(e); NA}))
}


#' Computing IC50 for a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param inflection : inflection point (EC50)
#' @param slope : Hill slope ( > 0 for decreasing curves)
#' @param minimum_dose : Minimum dose
#' @param maximum_dose : Maximum dose
#'
#' @return IC50 : The dose value where the curve intersects with y = 0.5, NA returned if they don't intersect in the given dose range.
#' @export
#'
#' @examples
compute_log_ic50 <- function(LL, UL, inflection, slope, minimum_dose, maximum_dose) {
  if((LL >= 0.5) | (UL <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (UL + (LL - UL)/(1 + (2^x/inflection)^slope)- 0.5)
    return(tryCatch(uniroot(f1, c(log2(minimum_dose), log2(maximum_dose)))$root,  error = function(x) NA))
  }
}

#' Computing mean square and median absolute errors for for a 4 parameter log-logistic dose response curve and the corresponding (dose,viability) pairs.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL : Upper asymptote
#' @param LL : Lower asymptote
#' @param slope : Hill slope
#' @param inflection : inflection point (EC50)
#'
#' @return List of mse and mad values.
#' @export
#'
#' @examples
compute_mse_mad <- function(FC, dose,  UL, LL,  slope, inflection) {
  FC.pred = UL  + (LL -UL )/(1 + (dose/inflection)^slope)
  residuals = FC - FC.pred
  return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
}


#' Fitting a dose response curve to the given dose and viability (FC) values.
#' This function fits 5 different dose-response functiosn to the given dose, viability pairs using dr4pl and drc packages and returns
#' the best one (lowest mse) among them.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL_low : Lower limit for the upper asymptote
#' @param UL_up : Upper limit for the upper asympotote
#' @param slope_decreasing: Should the curve to be constrained to be decreasing or not.
#'
#' @return Returns a single row data-frame with following columns:
#'          fit_name : Name of the fitting method with the highest explained variance (lowest mse)
#'          lower_limit : Lower asmpytote for the selected fit
#'          upper_limit : Upper asymptote for the selected fit
#'          slope : The Hill slope for the selected fit
#'          inflection : inflection point of the selected fit (EC50)
#'          mse : mse of the selected fit
#'          mad : mad of the selected fit
#'          frac_var_explained : Explained variance for the best fit
#'          successful_fit: If any of the fitting methods has positive frac_var_explained
#'          auc_riemann : The average measured viability value across doses
#'          minimum_dose : Minimum measured dose
#'          maximum_dose : Maximum measured dose
#'          auc : auc of the log-dose vs viability curve normalized to the measured dose-range (takes values between 0 and 1)
#'          log2_ic50 : Log2 of IC50 for the fitted curve
#' @export
#'
#' @examples
get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE, seed=42) {
  require(dr4pl)
  require(drc)
  require(tidyverse)
  require(magrittr)

  # Set seed for reproducibility of stochastic optimization
  set.seed(seed)


  # Fits a number of alternate models  to the DRC and chooses the best fit.

  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),
  #                          "lower_limit"=double(),
  #                          "upper_limit"=double(),
  #                          "slope"=double(),
  #                          "inflection"=double(),
  #                          "mse"=double(), "mad" =double(),
  #                          "frac_var_explained"=double())

  dose = as.numeric(dose)
  FC = as.numeric(FC)
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  riemann_auc <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
  var_data = var(FC)

  minimum_dose = min(dose); maximum_dose = max(dose)
  n_doses = length(unique(dose))
  results.df <- list(); ix = 1


  # FIT 1 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package


  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.

    results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                slope = -as.numeric(drc_model$coefficients[[1]]),
                                inflection = as.numeric(drc_model$coefficients[[4]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }



  # FIT 2 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package


  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    if((!slope_decreasing) | (as.numeric(drc_model$coefficients[[1]]) > 0)){
      mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.

      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                  upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                  slope = -as.numeric(drc_model$coefficients[[1]]),
                                  inflection = as.numeric(drc_model$coefficients[[4]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }


  # FIT 3 ---
  dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
                                                                       theta_3= -3, theta_4 = 0.01),
                                        lowerl = c(UL_low, -Inf, -Inf, 0),
                                        upperl = c(UL_up, Inf, slope_bound, 1.01),
                                        method.optim="Nelder-Mead"),
                                  error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )

  if (!dr4pl_initMan_optNM$convergence){
    if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) {
      dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
    }
  }

  if (dr4pl_initMan_optNM$convergence & all(is.finite(unlist(dr4pl_initMan_optNM$parameters)))){
    mse_mad <- compute_mse_mad(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                               dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])

    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                lower_limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }

  # FIT 4 ---
  dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})

  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
    }
  }


  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_mad <- compute_mse_mad(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                 dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                  lower_limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                  upper_limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                  slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                  inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }

  # FIT 5 ---
  dr4pl_initL <- tryCatch(dr4pl(dose, FC,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.01)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )

  if (!dr4pl_initL$convergence){
    if (!is.null(dr4pl_initL$dr4pl.robust)) {
      dr4pl_initL <- dr4pl_initL$dr4pl.robust
    }
  }

  if (dr4pl_initL$convergence & all(is.finite(unlist(dr4pl_initL$parameters)))){
    mse_mad <- compute_mse_mad(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                               dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])

    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                lower_limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }


  # Choose the best fit among the successful fits---
  results.df <- dplyr::bind_rows(results.df)

  if (nrow(results.df)>0){
    results.df <- results.df %>%
      dplyr::arrange(desc(frac_var_explained)) %>%
      head(1) %>%
      dplyr::mutate(successful_fit = TRUE,
                    auc_riemann = as.numeric(riemann_auc),
                    minimum_dose = minimum_dose, maximum_dose = maximum_dose) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(auc = compute_auc(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose),
                    log2_auc = log2(auc),
                    log2_ic50 = ifelse(n_doses>5, ## IC50 is not computed if we don't measure viability at at least 5 doses
                                       compute_log_ic50(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose),
                                       NA),
                    ic50= ifelse(is.na(log2_ic50), NA, 2^log2_ic50)) ## provide IC50 and log2IC50 only if we have enough doses and viability crosses0.5

  }else{
    results.df  <- data.frame(successful_fit=FALSE,
                              auc_riemann = as.numeric(riemann_auc),
                              minimum_dose = minimum_dose, maximum_dose = maximum_dose, 
                              auc = NA, log2_ic50 = NA, ic50 = NA)
  }

  return (results.df)
}




#' Creating the dose-response table for a single treatment viability screens.
#'
#' @param l2fc : The path to look for the log fold change input file.
#' @param out_path : Where to write the output file, make sure it ends with a "/". If left NULL it will write the output file to in_path.
#' @param cell_line_cols : List of columns to identify each cell line (model). Default: c("depmap_id", "cell_set", "pool_id").
#' @param treatment_cols : List of columns to idenfity each treatment. Default:  c("pert_id", "x_project_id", "pert_name", "pert_plate")
#' @param dose_column : Name of the column that represents the dose information. Default: "pert_dose"
#' @param l2fc_column : Name of the column that holds log2-viability values. Default: "l2fc"
#' @param type_column : Name of the column that contains the perturbation type. Only the rows that has value "trt_cp" will be kept in the calculations. Default: "pert_type"
#' @param cap_for_viability : The upper threshold for the viability values for before fitting the curves. Default is 1.5, any viability value above this value will be made equal to 1.5.
#'
#' @return It writes the calculated dose-response parameters to a file named as output_file_name in the out_path folder.
#' @export
#'
#' @examples create_drc_table(in_path = "data/sushi_io/testing_MTS-SEQ002-KF/test/")
create_drc_table <- function(LFC, screen_type,
                             cell_line_cols = c("depmap_id", "cell_set", "pool_id"),
                             treatment_cols = c("pert_id", "x_project_id", "pert_name", "pert_plate"),
                             dose_col = "pert_dose", l2fc_col = "median_l2fc", type_column = "pert_type",
                             cap_for_viability = 1.5, seed = 42) {
  require(data.table)
  require(tidyverse)
  require(rlang)

  # Check if the input table is populated
  if (nrow(LFC) == 0) {
    stop("LFC is empty!")
  }

  print(paste0("dose_column:", dose_col))
  print(paste0("l2fc_column:", l2fc_col))
  print(paste0("cell_line_cols:", cell_line_cols))
  print(paste0("treatment_cols:", treatment_cols))
  print(paste0("cap_for_viability:", cap_for_viability))
  print(paste0("Columns in LFC:", colnames(LFC)))

  # Check for missing critical columns
  necessary_columns <- unique(c(dose_col,l2fc_col,type_column,treatment_cols,cell_line_cols))

  print("Necessary columns:")
  print(necessary_columns)

  # Convert columns to numeric explicitly and check for NA generation
  LFC <- LFC %>%
    dplyr::mutate(
      dose_ = as.numeric(.data[[dose_col]]),
      l2fc_ = as.numeric(.data[[l2fc_col]])
    )

  # Verify conversion results
  print(paste0("NA in dose_: ", any(is.na(LFC$dose_))))
  print(paste0("NA in l2fc_: ", any(is.na(LFC$l2fc_))))

  # Ensure no NAs before proceeding
  if (any(is.na(LFC$dose_)) || any(is.na(LFC$l2fc_))) {
    stop("Conversion to numeric resulted in NA values!")
  }

  # Fit the curves for the grouped data
  if (screen_type == "CPS_SEQ") {
    DRC_SINGLE <- LFC %>%
      dplyr::filter(is.finite(dose_), is.finite(l2fc_)) %>%
      dplyr::group_by(across(all_of(c(cell_line_cols, treatment_cols)))) %>%
      dplyr::filter(length(unique(dose_)) > 4) %>%
      dplyr::summarise(if (dplyr::first(is_combination)) {
        # For combinations, use expanded fitting parameters
        get_best_fit(FC = pmin(2^l2fc_, cap_for_viability), dose = dose_,
                     UL_low = 0, UL_up = 1.01, slope_decreasing = FALSE, seed = seed)
      } else {
        # For non combinations, use default fitting parameters
        get_best_fit(FC = pmin(2^l2fc_, cap_for_viability), dose = dose_, seed = seed)
      }) %>%
      dplyr::ungroup()
  } else if (screen_type %in% c("MTS_SEQ", "EPS_SEQ", "APS_SEQ")) {
    DRC_SINGLE <- LFC %>%
      dplyr::filter(is.finite(dose_), is.finite(l2fc_)) %>%
      dplyr::group_by(across(all_of(c(cell_line_cols, treatment_cols)))) %>%
      dplyr::filter(length(unique(dose_)) > 4) %>%
      dplyr::summarise({
        get_best_fit(FC = pmin(2^l2fc_, cap_for_viability), dose = dose_, seed = seed)
      }) %>%
      dplyr::ungroup()
  } else {
    stop(screen_type, " is not a valid build type.")
  }
  return(DRC_SINGLE)
}
