library(useful)
library(tidyverse)
library(magrittr)
source("scripts/utils/kitchen_utensils.R")

# -----
# Load the raw data ----
# ----

filtered_counts <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqFilteredCounts.csv")


# ----
# Normalize relative to the spike-in control barcodes ----
# ----

source("scripts/normalize/normalize_functions.R")


CB_meta <- filtered_counts %>%
  dplyr::distinct(cb_ladder, cb_name, cb_log2_dose) %>%
  tidyr::drop_na()
  
normalized_counts = normalize(filtered_counts, 
                              id_cols = c("pcr_plate", "pcr_well"),
                              CB_meta = CB_meta, 
                              pseudocount = 0)

normalized_counts = add_pseudovalue(normalized_counts, c("pcr_plate", "pert_vehicle"),
                                    read_detection_limit = 10, negcon_type = "ctl_vehicle")

# -----
# QC : Note this section is OncRef specific?  ----
# ----

# FILTER NEGATIVE CONTROL WELLS THAT CONTROL BARCODES DOESN'T COVER AT LEAST 25% OF THE RANGE OF CELL LINE BARCODES OR
# AND ANY WELL THAT HAS SPEARMAN CORRELATION LESS THAN 0.85 WITH INTENDED DOSES. 

CB.quantiles <- normalized_counts %>% 
  dplyr::group_by(pcr_plate, pcr_well) %>% 
  dplyr::arrange(n) %>% 
  dplyr::mutate(quant = (1:n()) / n()) %>% 
  dplyr::filter(!is.na(cb_name)) %>% 
  dplyr::distinct(pcr_plate, pcr_well, pert_type, lua, n, quant) %>%
  dplyr::group_by(pcr_plate, pcr_well) %>% 
  dplyr::mutate(range = max(quant, na.rm = T) - min(quant, na.rm = T), n.b = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter((range >= 0.25) | (pert_type != "ctl_vehicle"), n.b > 5) 


# # Which wells are removed
# filtered_wells <- CB.quantiles %>% 
#   dplyr::left_join(normalized_counts) %>% 
#   dplyr::group_by(pcr_plate, pcr_well, pert_type) %>%
#   dplyr::summarise(concordance = cor(n, cb_log2_dose, use = "p", method = "spearman")) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::filter(concordance < 0.85)
# 
# 
# filtered_wells %>%  
#   dplyr::left_join(filtered_counts) %>% 
#   ggplot() +
#   geom_point(aes(x = log2(n), y = cb_log2_dose,
#                  color = pert_type)) + facet_wrap(paste0(pcr_plate, "/", pcr_well) ~ ., scales = "free")



normalized_counts <- CB.quantiles %>% 
  dplyr::left_join(normalized_counts) %>% 
  dplyr::group_by(pcr_plate, pcr_well, pert_type) %>%
  dplyr::summarise(concordance = cor(n, cb_log2_dose, use = "p", method = "spearman")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(concordance >= 0.85) %>% 
  dplyr::distinct(pcr_plate, pcr_well) %>%
  dplyr::left_join(normalized_counts)

rm(CB.quantiles) 



# DROP PLATES THAT DOESN'T HAVE AT LEAST 16 NEGCON WELLS
normalized_counts %<>% 
  dplyr::filter(pert_type == "ctl_vehicle") %>% 
  dplyr::distinct(pcr_plate, pcr_well) %>% 
  dplyr::count(pcr_plate) %>% 
  dplyr::filter(n > 16) %>% 
  dplyr::select(-n) %>% 
  dplyr::left_join(normalized_counts)



qc_table <- normalized_counts %>%
  dplyr::filter(!is.na(depmap_id)) %>% 
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon"), pool_id != "CTLBC") %>%
  dplyr::group_by(CompoundPlate, pcr_plate, cell_set, lua, pool_id, depmap_id) %>%
  dplyr::filter(is.finite(log2_normalized_n)) %>%  
  dplyr::summarise(
    median_negcon_raw_count = median(n[pert_type == "ctl_vehicle"], na.rm = T), 
    error_rate =  min(PRROC::roc.curve(scores.class0 = log2_normalized_n,
                                       weights.class0 = pert_type == "ctl_vehicle",
                                       curve = TRUE)$curve[,1] + 1 -
                        PRROC::roc.curve(scores.class0 = log2_normalized_n,
                                         weights.class0 = pert_type == "ctl_vehicle",
                                         curve = TRUE )$curve[,2])/2,
    NC.median = median(log2_normalized_n[pert_type == "ctl_vehicle"], na.rm = T),
    NC.mad = mad(log2_normalized_n[pert_type == "ctl_vehicle"], na.rm = T),
    PC.median = median(log2_normalized_n[pert_type == "trt_poscon"], na.rm = T),
    PC.mad = mad( log2_normalized_n[pert_type == "trt_poscon"], na.rm = T),
    NC.n = sum(pert_type == "ctl_vehicle", na.rm = T),
    PC.n = sum(pert_type == "trt_poscon", na.rm = T)) %>%
  dplyr::mutate(DR = NC.median - PC.median,
                SSMD = DR / sqrt(NC.mad ^ 2 + PC.mad ^ 2)) %>%
  dplyr::mutate(PASS = (error_rate <= 0.05) & (DR > 2) & (NC.n >= 16) & (PC.n >= 16) & (NC.mad <= 1) & (median_negcon_raw_count > 40)) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()




# ----
# Compute log2-fold-changes & correct for growth condition and negative control abundance.
# ----

source("scripts/compute_l2fc/compute_l2fc_functions.R")
source("scripts/bias_correction/bias_correction_functions.R")
source("scripts/collapse_replicates/collapse_replicates_functions.R")


l2fc = compute_l2fc(normalized_counts = normalized_counts,
                    control_type = "ctl_vehicle",
                    sig_cols = c("CompoundPlate", "SampleID","pert_dose", "pert_dose_unit","day"),
                    ctrl_cols = c("cell_set", "day", "pcr_plate", "replicate_plate"),
                    count_col_name = "log2_normalized_n",
                    cell_line_col = c("cell_set", "lua", "depmap_id", "pool_id", "growth_condition")) %>% 
  dplyr::semi_join(qc_table %>% 
                     dplyr::filter(PASS)) %>% 
  dplyr::mutate(negcon_log2_norm_n = log2(control_median_normalized_n)) %>%
  dplyr::group_split(bio_rep, CompoundPlate, SampleID, pert_dose, pert_dose_unit, day) %>%
  lapply(apply_bias_correction, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_condition") %>%
  dplyr::bind_rows() %>% 
  dplyr::ungroup()

# Collapse l2fcs ----

collapsed_l2fc = collapse_bio_reps(l2fc = l2fc, 
                                   sig_cols = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "day"), 
                                   cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"))



# -----
# Identify outlier pools and remove from l2fc tables  -----
# -----

dose_indices <- collapsed_l2fc %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose) %>% 
  dplyr::group_by(CompoundPlate, SampleID) %>% 
  dplyr::arrange(pert_dose) %>% 
  dplyr::mutate(dose_ix = 1:n()) %>%
  dplyr::ungroup()


monotonicity_flags <- dose_indices %>% 
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.prev = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix + 1)) %>%   
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix - 1)) %>%   
  dplyr::select(-dose_ix) %>% 
  dplyr::left_join(collapsed_l2fc %>%
                     dplyr::select(-median_l2fc_uncorrected, -num_bio_reps)) %>%  
  dplyr::left_join(collapsed_l2fc %>%
                     dplyr::select(-median_l2fc_uncorrected, -num_bio_reps) %>% 
                     dplyr::rename(pert_dose.prev = pert_dose,  median_l2fc.prev = median_l2fc)) %>%  
  dplyr::left_join(collapsed_l2fc %>%                
                     dplyr::select(-median_l2fc_uncorrected, -num_bio_reps) %>% 
                     dplyr::rename(pert_dose.next = pert_dose,  median_l2fc.next = median_l2fc)) %>%  
  dplyr::mutate(flag.down = median_l2fc < -2, flag.up =  median_l2fc > -1,
                flag.down.prev = ifelse(is.na(pert_dose.prev), FALSE, median_l2fc.prev < -2), 
                flag.up.prev =  ifelse(is.na(pert_dose.prev), TRUE, median_l2fc.prev > -1),
                flag.down.next = ifelse(is.na(pert_dose.next), TRUE, median_l2fc.next < -2),
                flag.up.next =  ifelse(is.na(pert_dose.next), FALSE, median_l2fc.next > -1)) %>% 
  dplyr::mutate(flag1 = flag.down & flag.up.next & flag.up.prev,
                flag2 = flag.up & flag.down.next & flag.down.prev)

outlier.trt.pools <- monotonicity_flags %>% 
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, cell_set, pool_id) %>% 
  dplyr::summarise(n.f1 = mean(flag1, na.rm = T),
                   n.f2 = mean(flag2, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(pmax(n.f1, n.f2) > 0.25) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, cell_set, pool_id) 

if(nrow(outlier.trt.pools) > 0){
  l2fc %<>% 
    dplyr::anti_join(outlier.trt.pools)
  
  collapsed_l2fc %<>% 
    dplyr::anti_join(outlier.trt.pools)
}

# -----
# Fitting dose-response curves ----
# -----

source("scripts/drc/dose_response_functions.R")

# drc_table <- l2fc %>% 
#   dplyr::group_split(CompoundPlate, SampleID) %>% 
#   lapply(create_drc_table,
#          l2fc_col = "l2fc",
#          screen_type = "MTS_SEQ", cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
#          treatment_cols = c("CompoundPlate", "SampleID"), dose_col = "pert_dose", type_col = "pert_type", cap_for_viability = 1.5) %>% 
#   dplyr::bind_rows() %>% 
#   dplyr::filter(successful_fit)
#
#
# drc_table_uncorrected <- l2fc %>% 
#   dplyr::group_split(CompoundPlate, SampleID) %>% 
#   lapply(create_drc_table,
#          l2fc_col = "l2fc_uncorrected",
#          screen_type = "MTS_SEQ", cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
#          treatment_cols = c("CompoundPlate", "SampleID"), dose_col = "pert_dose", type_col = "pert_type", cap_for_viability = 1.5) %>% 
#   dplyr::bind_rows() %>% 
#   dplyr::filter(successful_fit)
# 
# drc_table <- drc_table %>% 
#   dplyr::mutate(response = "corrected") %>% 
#   dplyr::bind_rows(drc_table_uncorrected %>% 
#                      dplyr::mutate(response = "uncorrected"))


plates = l2fc$CompoundPlate %>% unique()
drc_table <- list(); drc_table_uncorrected <- list(); 
ix = 1;   now = Sys.time()
for(plate in plates){
  print(plate)
  print(Sys.time() - now)
  now = Sys.time()
  
  drc_table[[ix]] <- l2fc %>% 
    dplyr::filter(CompoundPlate == plate) %>% 
    dplyr::group_split(SampleID) %>% 
    parallel::mclapply(create_drc_table,
                       l2fc_col = "l2fc",
                       screen_type = "MTS_SEQ", cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
                       treatment_cols = c("CompoundPlate", "SampleID"), dose_col = "pert_dose", type_col = "pert_type", cap_for_viability = 1.5,
                       mc.cores = pmax(parallel::detectCores() - 2,1)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(successful_fit)
  
  print(plate)
  print(Sys.time() - now)
  now = Sys.time()
  
  drc_table_uncorrected[[ix]] <- l2fc %>% 
    dplyr::filter(CompoundPlate == plate) %>% 
    dplyr::group_split(SampleID) %>% 
    parallel::mclapply(create_drc_table,
           l2fc_col = "l2fc_uncorrected",
           screen_type = "MTS_SEQ", cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
           treatment_cols = c("CompoundPlate", "SampleID"), dose_col = "pert_dose", type_col = "pert_type", cap_for_viability = 1.5,
           mc.cores = pmax(parallel::detectCores() - 2,1)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(successful_fit)
  

  ix <- ix + 1 
}

drc_table <- drc_table %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(response = "corrected") %>% 
  dplyr::bind_rows(drc_table_uncorrected %>% 
                     dplyr::bind_rows() %>% 
                     dplyr::mutate(response = "uncorrected")) 



drc_table %>% 
  write_csv("OncRef data/drc_table.csv")



# Compute fitted l2fc values ----

fitted_l2fc <- collapsed_l2fc %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, pert_dose_unit) %>% 
  dplyr::inner_join(drc_table) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(l2fc_fitted = log2(lower_limit + (upper_limit - lower_limit) / (1 + (pert_dose / inflection)^-slope))) %>%  
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, pert_dose_unit, cell_set, depmap_id, pool_id, lua, l2fc_fitted, response) %>% 
  dplyr::mutate(response = ifelse(response == "corrected", "l2fc_fitted", "l2fc_uncorrected_fitted")) %>% 
  tidyr::pivot_wider(names_from = response, values_from = l2fc_fitted) 
  
  
collapsed_l2fc %<>% 
  dplyr::left_join(fitted_l2fc) 


# SAVE FILES ----

qc_table %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqQC.csv")

l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFC.csv")

collapsed_l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFCCollapsed.csv")

drc_table %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqDRC.csv")


# Portal files 


compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundList.csv")


priority_table <- qc_table %>%
  dplyr::filter(PASS) %>% 
  dplyr::group_by(depmap_id, pool_id, lua, cell_set, CompoundPlate) %>% 
  dplyr::summarize(NC.mad = median(NC.mad, na.rm = T)) %>% 
  dplyr::ungroup()

priority_table <- drc_table %>% 
  dplyr::distinct(CompoundPlate, SampleID, depmap_id, pool_id, lua, cell_set) %>% 
  dplyr::left_join(priority_table) %>% 
  dplyr::group_by(CompoundPlate, SampleID, depmap_id) %>% 
  dplyr::arrange(NC.mad) %>% dplyr::mutate(rank = 1:n()) %>% dplyr::filter(rank == 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-rank, -NC.mad) %>% 
  dplyr::distinct() 


drc_table %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(Prioritized, successful_fit) %>% 
  dplyr::inner_join(priority_table) %>% 
  dplyr::rename(ModelID = depmap_id, EC50 = inflection, LowerAsymptote = lower_limit, UpperAsymptote = upper_limit, Slope = slope) %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqResponseCurves.csv")



drc_table %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter( Prioritized, successful_fit) %>% 
  dplyr::inner_join(priority_table) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2_auc") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqLog2AUCMatrix.csv")


drc_table %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter( Prioritized, successful_fit) %>% 
  dplyr::inner_join(priority_table) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "auc") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqAUCMatrix.csv")


drc_table %>% 
  dplyr::filter(response == "uncorrected") %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter( Prioritized, successful_fit) %>% 
  dplyr::inner_join(priority_table) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2_auc") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqLog2AUCMatrixUncorrected.csv")


drc_table %>% 
  dplyr::filter(response == "uncorrected") %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter( Prioritized, successful_fit) %>% 
  dplyr::inner_join(priority_table) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "auc") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqAUCMatrixUncorrected.csv")



collapsed_l2fc %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(Prioritized) %>% 
  dplyr::inner_join(priority_table) %>%  
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, depmap_id, l2fc_fitted) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "l2fc_fitted") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqLog2ViabilityCollapsedMatrix.csv")


collapsed_l2fc %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(Prioritized) %>% 
  dplyr::inner_join(priority_table) %>%  
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit) %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLog2ViabilityCollapsedConditions.csv")



collapsed_l2fc %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(Prioritized) %>% 
  dplyr::inner_join(priority_table) %>%  
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_uncorrected_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, depmap_id, l2fc_uncorrected_fitted) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "l2fc_uncorrected_fitted") %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqLog2ViabilityCollapsedMatrixUncorrected.csv")


collapsed_l2fc %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(Prioritized) %>% 
  dplyr::inner_join(priority_table) %>%  
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_uncorrected_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit) %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLog2ViabilityCollapsedConditionsUncorrected.csv")




LFC_ <- compound_annotations %>% 
  dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized) %>% 
  dplyr::filter(Prioritized) %>% 
  dplyr::left_join(l2fc) %>% 
  dplyr::left_join(qc_table) %>% 
  dplyr::inner_join(priority_table) %>% 
  dplyr::filter(is.finite(l2fc), PASS) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate, l2fc) %>% 
  dplyr::group_by(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::mutate(Replicate = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) 





LFC.conditions <- LFC_ %>% 
  dplyr::distinct(SampleID, Dose, DoseUnit, CompoundPlate, Replicate) %>% 
  dplyr::arrange(SampleID, CompoundPlate, Dose, DoseUnit, Replicate) %>% 
  dplyr::mutate(Label = 1:n() - 1) 



LFC_ %>% 
  dplyr::left_join(LFC.conditions) %>% 
  dplyr::mutate(viability = pmin(2^l2fc, 1.5)) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "viability") %>%
  write.csv("OncRef data/PRISMOncologyReferenceSeqViabilityMatrix.csv")

LFC.conditions %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqViabilityConditions.csv")






conf.pools <- filtered_counts %>% 
  dplyr::distinct(pool_id, depmap_id, cell_set) %>% 
  dplyr::filter(!is.na(depmap_id)) %>% 
  dplyr::mutate(ix = 1) %>%
  reshape2::acast(depmap_id ~ pool_id + cell_set, value.var = "ix", fill = 0) 


conf.pools.2 <- filtered_counts %>% 
  dplyr::distinct(depmap_id, cell_set) %>% 
  dplyr::filter(!is.na(depmap_id)) %>% 
  dplyr::mutate(ix = 1) %>%
  reshape2::acast(depmap_id ~ cell_set, value.var = "ix", fill = 0) 

conf.pools <- cbind(conf.pools, conf.pools.2[rownames(conf.pools), ])

conf.pools <- conf.pools[, !duplicated(t(conf.pools))]


conf.QC <- qc_table %>% 
  dplyr::inner_join(priority_table) %>% 
  dplyr::distinct(depmap_id, NC.median, NC.mad, PC.median, PC.mad, DR, SSMD) %>% 
  tidyr::pivot_longer(cols = c(2:7)) %>% 
  dplyr::filter(is.finite(value)) %>% 
  reshape2::acast(depmap_id ~ name, fun.aggregate = median)


conf.pools %>% 
  reshape2::melt() %>% 
  dplyr::bind_rows(conf.QC %>% 
                     reshape2::melt()) %>% 
  reshape2::acast(Var1 ~ Var2) %>% 
  write.csv("OncRef data/PRISMOncologyReferenceSeqConfounderMatrix.csv")
