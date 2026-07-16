library(tidyverse)
library(parallel)
library(matrixStats) # Added for fast matrix operations

source("scripts/oncref_data_processing/manuscript_biomarker_utilities.r")
file_path <- "~/code/depmap_data/depmap_26Q1_internal.h5" # !!!

set.seed(23)


# -----
# 1. Load the PRISM data ----
# -----

CompoundList <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceCompoundList.csv") 


LAUC.harmonized <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceHarmonizedLog2AUCMatrix.csv") %>% 
  tibble::column_to_rownames("V1") %>% 
  as.matrix() 




# -----
# 2. Target-recovery functions -----
# -----

bm.lauc <- target_recovery(Y = LAUC.harmonized, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.lauc, "OncRef data/PRISMOncologyReferenceHarmonizedLog2AUCUnivariateBiomarkers.csv")

bm.lauc %>% 
  dplyr::group_by(feature_set) %>% 
  dplyr::filter(is.target) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(CompoundName, Readout, feature_set) %>%
  dplyr::mutate(mr = min(rank[is.target])) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(CompoundName, Readout, feature_set, mr) %>% 
  dplyr::mutate(mr = pmin(mr, 100)) %>% 
  ggplot() +
  stat_ecdf(aes(x = mr, color = Readout)) +
  scale_x_log10() +
  facet_wrap(feature_set ~ .)
  
  
# ----
# 3. Load DepMap data into a single matrix (Optimized) ----
# ----

datasets_to_load <- c(
  "CRISPR" = "CRISPR", 
  "RNAi" = "RNAi", 
  "Expression" = "EXP", 
  "CopyNumber" = "CN", 
  "Mutation" = "MUT", 
  "Lineage" = "LIN", 
  "Fusion" = "FUS"
)

target_rows <- rownames(LAUC.harmonized)

# Loop over datasets, clean, align to LAUC, and impute medians natively as matrices
matrix_list <- lapply(names(datasets_to_load), function(ds_name) {
  prefix <- datasets_to_load[[ds_name]]
  mat <- read_dataset(file = file_path, dataset = ds_name)
  
  # Intersect with LAUC initially
  cl <- intersect(rownames(mat), target_rows)
  mat <- mat[cl, , drop = FALSE]
  
  # Keep columns with > 90% finite values
  mat <- mat[, colMeans(is.finite(mat)) > 0.9, drop = FALSE]
  colnames(mat) <- paste0(prefix, "_", colnames(mat))
  
  # Initialize an empty matrix perfectly aligned to LAUC rows
  aligned_mat <- matrix(NA_real_, nrow = length(target_rows), ncol = ncol(mat), 
                        dimnames = list(target_rows, colnames(mat)))
  
  # Fill existing rows
  common_rows <- intersect(target_rows, rownames(mat))
  aligned_mat[common_rows, ] <- mat[common_rows, ]
  
  # Fast median imputation for missing values
  col_meds <- matrixStats::colMedians(aligned_mat, na.rm = TRUE)
  na_idx <- which(is.na(aligned_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    aligned_mat[na_idx] <- col_meds[na_idx[, 2]]
  }
  
  return(aligned_mat)
})

# Column-bind all matrices instantly (Replaces reshape2::melt -> join -> acast)
X <- do.call(cbind, matrix_list)

# Filter final matrix by variance
X <- X[, matrixStats::colVars(X, na.rm = TRUE) > 0.005, drop = FALSE]


# -----
# 4. RANDOM FOREST MODELS -----
# -----

RF.lauc <- biomarker_suite_rf_cv(
  X, LAUC.harmonized, 
  biomarker_file = file_path, 
  CompoundList = CompoundList,                                       
  bm_th = 0.05, bm_R = 10, bm_R2 = 50,
  K = 10,
  seed = 23
)

saveRDS(RF.lauc, "OncRef data/PRISMOncologyReferenceHarmonizedLog2AUCMultivariateBiomarkers.RDS") 


# ----
# 5. BIOMARKER SUMMARY TABLES ----
# ----

# Auxiliary tables
BM.Summary.Table <- RF.lauc$model_performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(model, cn, CompoundName) %>%  
  dplyr::summarise(mse = mean(mse), r.sd = sd(r), r = mean(r), var.y = mean(var.y.test), .groups = "drop_last") %>% 
  dplyr::mutate(r2 = 1 - mse / var.y) %>% 
  dplyr::group_by(cn) %>%
  dplyr::mutate(
    r.m = max(r[!model %in% c("targets", "extended")], na.rm = TRUE),
    n.t = length(setdiff(model, c("targets", "extended"))),
    model.class = dplyr::case_when(
      model == "extended" ~ "Extended",
      model == "targets" ~ "Targets",
      r == r.m ~ "Best Single Target",
      TRUE ~ "Other Targets"
    )
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-r.m)




# -----



DF <- RF.lauc$predictions %>% 
  dplyr::filter(type == "test") %>%
  dplyr::select(-y.hat.oob) %>% 
  dplyr::distinct() %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(CompoundName, model) %>% 
  dplyr::arrange(y.hat) %>% 
  dplyr::mutate(
    n = dplyr::row_number(), 
    N = dplyr::n(),
    p = var(y) * (1 / n + 1 / (N - n)), 
    cs = cumsum(y), 
    s = sum(y),
    m1 = cs / n, 
    m2 = (s - cs) / (N - n),
    t = -(m1 - m2) / sqrt(p)
  ) %>%
  dplyr::select(-p, -s, -cs, -m1, -m2, -K) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct()


DF <- BM.Summary.Table %>% 
  dplyr::left_join(
    DF %>% 
      dplyr::filter(is.finite(t)) %>% 
      dplyr::group_by(CompoundName, model, N) %>% 
      dplyr::summarize(
        t.mean = mean(t),
        t.peak = max(t),
        n.peak = min(n[t == t.peak]),
        .groups = "drop"
      ),
    by = c("CompoundName", "model")
  ) %>% 
  dplyr::left_join(
    CompoundList %>% dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism),
    by = "CompoundName"
  )


# Scores table - Pivot directly without splitting the dataframe
Scores.Table <- BM.Summary.Table %>% 
  dplyr::filter(model.class != "Other Targets") %>% 
  dplyr::distinct(CompoundName, cn, n.t, r, r.sd, model.class) %>%
  tidyr::pivot_wider(
    names_from = model.class, 
    values_from = c(r, r.sd),
    names_glue = "{.value}_{model.class}" # Formats as r_Targets, r.sd_Targets, etc.
  ) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    PolypharmacologyScore = ifelse(n.t > 1, (r_Targets - `r_Best Single Target`) / sqrt((`r.sd_Targets`^2 + `r.sd_Best Single Target`^2) / 10), 0),
    ExcessPredictabilityScore = ifelse(PolypharmacologyScore > 0, 
                                       (r_Extended - r_Targets) / sqrt((`r.sd_Targets`^2 + `r.sd_Extended`^2) / 10),
                                       (r_Extended - `r_Best Single Target`) / sqrt((`r.sd_Extended`^2 + `r.sd_Best Single Target`^2) / 10))
  ) %>% 
  dplyr::mutate(
    PolypharmacologyScore = pmax(PolypharmacologyScore, 0),
    ExcessPredictabilityScore = pmax(ExcessPredictabilityScore, 0),
    Best.r = pmax(r_Extended, r_Targets, `r_Best Single Target`, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(CompoundName, cn, n.t, PolypharmacologyScore, ExcessPredictabilityScore, Best.r, r_Extended, r_Targets, `r_Best Single Target`)


Scores.Table <- DF %>% 
  dplyr::filter(model %in% c("targets", "extended")) %>% 
  dplyr::distinct(cn, CompoundName, model, t.mean, t.peak, N, n.peak) %>% 
  dplyr::mutate(SelectivityScore = t.peak - pmax(t.mean, 0)) %>% 
  tidyr::pivot_wider(names_from = model, values_from = c(SelectivityScore, t.mean, t.peak, N, n.peak)) %>% 
  dplyr::left_join(Scores.Table, by = c("cn", "CompoundName")) %>% 
  dplyr::rename(
    OnTargetPolypharmacologyScore = PolypharmacologyScore,
    OffTargetPolypharmacologyScore = ExcessPredictabilityScore,
    n.targets = n.t
  ) %>% 
  dplyr::select(cn, CompoundName, Best.r,
                OnTargetPolypharmacologyScore, OffTargetPolypharmacologyScore,
                SelectivityScore_extended, SelectivityScore_targets,
                r_Extended, r_Targets, `r_Best Single Target`,
                n.peak_extended, t.peak_extended, t.mean_extended, N_extended, 
                n.peak_targets, t.peak_targets, t.mean_targets, N_targets,
                n.targets)


# Variable importances 
Importance.Table <- RF.lauc$variable_importances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(cn, model, K) %>% 
  dplyr::mutate(imp = imp / sum(imp)) %>%  
  dplyr::group_by(cn, CompoundName, model, var) %>% 
  dplyr::summarise(imp = sum(imp) / 10, .groups = "drop_last") %>% 
  dplyr::arrange(desc(imp)) %>% 
  dplyr::mutate(rank = dplyr::row_number()) %>% 
  dplyr::ungroup() 

# Predictability results
Predictability.Table <- RF.lauc$model_performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(cn, CompoundName, model) %>% 
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  dplyr::select(cn, CompoundName, model, mse, r2, r, var.y.test)


# ----
# 6. SAVE RESULTS ----
# ----

readr::write_csv(Predictability.Table, "OncRef data/PRISMOncologyReferenceHarmonizedLog2AUCMultivariateBiomarkersModelTable.csv")
readr::write_csv(Importance.Table, "OncRef data/PRISMOncologyReferenceHarmonizedLog2AUCMultivariateBiomarkersVarImpTable.csv")
readr::write_csv(Scores.Table, "OncRef data/PRISMOncologyReferenceHarmonizedLog2AUCMultivariateBiomarkersScoresTable.csv")





# ----
# TO-DO: Repeat for transformed Luminex, LFC and LAUC.

compound_annotations <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceCompoundList.csv")  %>% 
  dplyr::filter(Prioritized, hasData,
                Readout == "Luminex") %>% 
  dplyr::ungroup() 


drc_table_lum_harmonized <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceHarmonizedLumDRC.csv") 


transformed_lum_LAUC <- drc_table_lum_harmonized %>% 
  dplyr::inner_join(compound_annotations) %>% 
  dplyr::filter(successful_fit) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2_auc")
  
  
bm.transformed_lauc <- target_recovery(Y = transformed_lum_LAUC, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.transformed_lauc, "OncRef data/PRISMOncologyReferenceHarmonizedLumLog2AUCUnivariateBiomarkers.csv")


collapsed_l2fc_lum_harmonized <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceHarmonizedLumLFCCollapsed.csv")

transformed_lfc <- collapsed_l2fc_lum_harmonized %>%
  dplyr::inner_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate)) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(SampleID, "::", pert_dose)) %>% 
  dplyr::distinct(Label, SampleID, depmap_id, l2fc_fitted) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "l2fc_fitted") 


# Note these are fitted ones 
bm.transformed_lfc <- target_recovery(Y = transformed_lfc, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.transformed_lfc, "OncRef data/PRISMOncologyReferenceHarmonizedLumLFCUnivariateBiomarkers.csv")




# -----
lum_lfc <- collapsed_l2fc_lum_harmonized %>%
  dplyr::inner_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate)) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, median_l2fc, CompoundName) %>% 
  dplyr::mutate(Label = paste0(SampleID, "::", pert_dose)) %>% 
  dplyr::distinct(Label, SampleID, depmap_id, median_l2fc) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "median_l2fc") 



bm.lum_lfc <- target_recovery(Y = lum_lfc, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.lum_lfc, "OncRef data/PRISMOncologyReferenceLumLFCUnivariateBiomarkers_notfitted.csv")



lum.tr_lfc <- collapsed_l2fc_lum_harmonized %>%
  dplyr::inner_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate)) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, median_l2fc_transformed, CompoundName) %>% 
  dplyr::mutate(Label = paste0(SampleID, "::", pert_dose)) %>% 
  dplyr::distinct(Label, SampleID, depmap_id, median_l2fc_transformed) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "median_l2fc_transformed") 



bm.lum.tr_lfc <- target_recovery(Y = lum.tr_lfc, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.lum.tr_lfc, "OncRef data/PRISMOncologyReferenceHarmonizedLumLFCUnivariateBiomarkers_notfitted.csv")



# Comparison ----

comp_table <- bm.lum.tr_lfc %>% 
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, feature, feature_set, CompoundName, rank, correlation_coef) %>% 
  dplyr::full_join(bm.lum_lfc %>% 
                     dplyr::filter(is.target) %>% 
                     dplyr::distinct(cn, feature, feature_set, CompoundName, rank, correlation_coef),
                   by = c("cn", "feature", "feature_set", "CompoundName")) 



comp_table %<>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(rank.x = ifelse(is.finite(rank.x), rank.x, 100),
                rank.y = ifelse(is.finite(rank.y), rank.y, 100)) %>% 
  dplyr::ungroup()
  
  

comp_table %>%
  dplyr::group_by(CompoundName, feature, feature_set) %>% 
  dplyr::summarise(rank.x = min(rank.x), rank.y = min(rank.y)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(abs(rank.x - rank.y) > 10, abs(log2(rank.x / rank.y)) > 1,
                pmin(rank.x, rank.y) < 10) %>%
  dplyr::mutate(rr = ifelse(rank.x < rank.y, "P", "N")) %>%
  View
  dplyr::count(feature_set, rr)  %>% 
  tidyr::pivot_wider(names_from = rr, values_from = n, values_fill = 0)

comp_table %>% 
    dplyr::filter(pmin(rank.x, rank.y) < 10) %>%
    dplyr::mutate(rr = ifelse(rank.x < rank.y, "P", "N")) %>% 
    dplyr::count(feature_set, rr)  %>% 
    tidyr::pivot_wider(names_from = rr, values_from = n, values_fill = 0) 

comp_table %>% 
  dplyr::filter(abs(rank.x - rank.y) > 10, #abs(log2(rank.x / rank.y)) > 1,
                pmin(rank.x, rank.y) < 10) %>%
  dplyr::mutate(rr = ifelse(rank.x < rank.y, "P", "N")) %>% 
  dplyr::count(feature_set, rr)  %>% 
  tidyr::pivot_wider(names_from = rr, values_from = n, values_fill = 0)
  


comp_table %>% 
  dplyr::group_by(CompoundName, feature, feature_set) %>% 
  dplyr::summarise(rank.x = min(rank.x), rank.y = min(rank.y),
                   r2.x= max(correlation_coef.x^2, na.rm = T),  
                   r2.y= max(correlation_coef.y^2, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(abs(rank.x - rank.y) > 10, #abs(log2(rank.x / rank.y)) > 1,
                pmin(rank.x, rank.y) < 10) %>%
  dplyr::mutate(rr = ifelse(rank.x < rank.y, "P", "N")) %>% 
  ggplot() + 
  geom_point(aes(x = r2.x, y = r2.y,
                 color = rr)) +
  geom_abline() + 
  facet_wrap(feature_set ~ .)
  




bm.lum_lfc %>%  
  dplyr::filter(rank <= 1) %>% 
  dplyr::count(feature, feature_set, status) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::rename(n.orig = n) %>% 
  dplyr::full_join(bm.lum.tr_lfc %>%  
                     dplyr::filter(rank <= 1) %>% 
                     dplyr::count(feature, feature_set, status) %>% 
                     dplyr::arrange(desc(n))) %>%
  dplyr::mutate(n.orig = ifelse(is.finite(n.orig), n.orig, 0),
                n = ifelse(is.finite(n), n, 0)) %>% 
  dplyr::filter(n.orig > 10 | n > 10,
                feature_set != "Lineage") %>% 
  View
  ggplot() +
  geom_point(aes(x = 1 + n.orig, y = 1 + n, color = feature_set)) +
  geom_abline() + 
  facet_wrap(feature_set ~ .)
  scale_x_log10() + scale_y_log10()
    