library(useful)
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
source("scripts/utils/kitchen_utensils.R")

# -----
# Load the raw data ----
# ----

filtered_counts <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqFilteredCounts.csv") %>% 
  dplyr::bind_rows(data.table::fread("~/Downloads/all_ADC_AIR_Seq_filtered_counts_sequencing.csv")) %>% 
  dplyr::bind_rows(data.table::fread("OncRef data/air/air001_filtered_counts.csv") %>% 
                     dplyr::rename(CompoundPlate = pert_plate, SampleID = pert_id)) %>% 
  dplyr::distinct() # %>% 
  # dplyr::mutate(growth_condition = ifelse(growth_condition == "Suspension", "S", "A")) 
  
# red <- data.table::fread("OncRef data/air/redactions.csv")

# filtered_counts %<>% 
#   dplyr::anti_join(red)

data.table::fread("~/Downloads/List of antibodies, ADCs, & payloads screened in PRISM - AIR agents for DMC 26Q3 release.csv") %>% 
  dplyr::filter(SampleID != "") %>% 
  dplyr::left_join(filtered_counts %>%
                     dplyr::distinct(CompoundPlate, SampleID) %>% dplyr::mutate(mustafa = T)) %>%
  dplyr::mutate(mustafa = !is.na(mustafa)) %>%
  dplyr::count(Screen, CompoundPlate, mustafa)
  View

filtered_counts %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_type) %>% 
  dplyr::filter(CompoundPlate %in% c("PAPS006", "PAPS007", "PAPS008")) %>% 
  dplyr::full_join(data.table::fread("~/Downloads/List of antibodies, ADCs, & payloads screened in PRISM - AIR agents for DMC 26Q3 release.csv") %>% 
                     dplyr::filter(SampleID != "") %>% 
                     dplyr::distinct(CompoundPlate, SampleID, Screen)) %>% 
  View
  
  dplyr::count(CompoundPlate)
  head

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
# QC : Note this section is OncRef specific!  ----
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
# 


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
  dplyr::group_by(CompoundPlate, cell_set, lua, pool_id, depmap_id) %>%
  dplyr::mutate(n.PASS = sum(PASS, na.rm = T)) %>% 
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
                     dplyr::filter(PASS, n.PASS > 1)) %>% 
  dplyr::mutate(negcon_log2_norm_n = log2(control_median_normalized_n)) %>%
  dplyr::group_split(bio_rep, CompoundPlate, SampleID, pert_dose, pert_dose_unit, day) %>%
  lapply(apply_bias_correction, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_condition") %>%
  dplyr::bind_rows() %>% 
  dplyr::ungroup() %>% 
  apply_shrinkage_correction_mild(shrunk_col_name = "l2fc_adaptive_mild", batch_vars = c("growth_condition") ) %>% 
  dplyr::rename(lambda_mild = lambda) %>% 
  apply_shrinkage_correction_heavy(shrunk_col_name = "l2fc_adaptive_heavy", batch_vars = c("growth_condition") ) %>% 
  dplyr::rename(lambda_heavy = lambda) %>% 
  apply_shrinkage_correction_mild(shrunk_col_name = "l2fc_adaptive_mild_complete") %>% 
  dplyr::rename(lambda_mild_complete = lambda) %>% 
  apply_shrinkage_correction_heavy(shrunk_col_name = "l2fc_adaptive_heavy_complete") %>% 
  dplyr::rename(lambda_heavy_complete = lambda) 


# Collapse l2fcs ----

collapsed_l2fc = collapse_bio_reps(l2fc = l2fc, 
                                   median_cols = c("l2fc", "l2fc_uncorrected", "l2fc_adaptive_mild", "l2fc_adaptive_heavy",  "l2fc_adaptive_heavy_complete",  "l2fc_adaptive_mild_complete"),
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
                     dplyr::select(-median_l2fc_uncorrected, -median_l2fc_adaptive_mild,-median_l2fc_adaptive_heavy, -median_l2fc_adaptive_mild_complete, -median_l2fc_adaptive_heavy_complete, -num_bio_reps)) %>%  
  dplyr::left_join(collapsed_l2fc %>%
                     dplyr::select(-median_l2fc_uncorrected, -median_l2fc_adaptive_mild,-median_l2fc_adaptive_heavy, -median_l2fc_adaptive_mild_complete, -median_l2fc_adaptive_heavy_complete, -num_bio_reps) %>% 
                     dplyr::rename(pert_dose.prev = pert_dose,  median_l2fc.prev = median_l2fc)) %>%  
  dplyr::left_join(collapsed_l2fc %>%                
                     dplyr::select(-median_l2fc_uncorrected, -median_l2fc_adaptive_mild, -median_l2fc_adaptive_heavy, -median_l2fc_adaptive_mild_complete, -median_l2fc_adaptive_heavy_complete, -num_bio_reps) %>% 
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


collapsed_l2fc %>%
  dplyr::inner_join(outlier.trt.pools %>%
                      dplyr::rename(pert_dose2 = pert_dose)) %>%
  ggplot(aes(x = log2(pert_dose), y = median_l2fc, group = depmap_id,
             )) +
  geom_line(lwd = 0.25) +
  geom_point(aes(color = (pert_dose2 == pert_dose)), show.legend = F) +
  facet_wrap(pool_id ~ CompoundPlate + SampleID)



if(nrow(outlier.trt.pools) > 0){
  l2fc %<>% 
    dplyr::anti_join(outlier.trt.pools)
  
  collapsed_l2fc %<>% 
    dplyr::anti_join(outlier.trt.pools)
}

# -----
# SAVE FILES ----
# -----


l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFC_adaptive_all.csv")

collapsed_l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFCCollapsed_adaptive_all.csv")



# Portal files !!!! ------


compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundList.csv") %>%
  dplyr::bind_rows(data.table::fread("~/Downloads/Ab_ADCs for manuscript - Detailed metadata.csv") %>% 
                     dplyr::distinct(CompoundPlate, CompoundName, SampleID, TargetOrMechanism, GeneSymbolOfTargets, Synonyms) %>% 
                     dplyr::mutate(Prioritized = TRUE) ) %>%
  dplyr::bind_rows(data.table::fread("~/Downloads/List of antibodies, ADCs, & payloads screened in PRISM - All Ab & ADCs screened to date (as of 3_2026).csv") %>% 
                     dplyr::select(2, 5,6, 9) %>% 
                     dplyr::rename(CompoundName = short_name, CompoundPlate = plate, SampleID = compound_id, GeneSymbolOfTargets = ab_target_gene)) %>% 
  dplyr::distinct(CompoundName, CompoundPlate, SampleID, GeneSymbolOfTargets)




compound_annotations %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqCompoundListExtended.csv")



lfc_df <- collapsed_l2fc %>%
  dplyr::inner_join(compound_annotations) 



file = "~/code/depmap_data/depmap_26Q1_internal.h5"
source("scripts/biomarker/biomarker_functions.R")


lfc_df %<>% 
  dplyr::mutate(cn = paste0(CompoundPlate, "::", SampleID, "::", pert_dose))
  
  
lfc.m.c <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc", fun.aggregate = median)


lfc.m.u <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_uncorrected)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_uncorrected", fun.aggregate = median)
  

lfc.m.am <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_mild)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_mild", fun.aggregate = median)


lfc.m.ah <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_heavy)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_heavy", fun.aggregate = median)


lfc.m.amc <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_mild)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_mild_complete", fun.aggregate = median)


lfc.m.ahc <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_heavy)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_heavy_complete", fun.aggregate = median)




bm.corrected <- univariate_biomarker_table(lfc.m.c, file = file,  v.X.min = 0.005,
                           features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                           regression_coef = F, stability_score = F, rank.max = 100)

bm.uncorrected <- univariate_biomarker_table(lfc.m.u, file = file, v.X.min = 0.005,
                                           features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                                           regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.m <- univariate_biomarker_table(lfc.m.am, file = file, v.X.min = 0.005,
                                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                                             regression_coef = F, stability_score = F, rank.max = 100)

bm.adaptive.h <- univariate_biomarker_table(lfc.m.ah, file = file, v.X.min = 0.005,
                                            features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                                            regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.mc <- univariate_biomarker_table(lfc.m.amc, file = file, v.X.min = 0.005,
                                            features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                                            regression_coef = F, stability_score = F, rank.max = 100)

bm.adaptive.hc <- univariate_biomarker_table(lfc.m.ahc, file = file, v.X.min = 0.005,
                                            features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                                            regression_coef = F, stability_score = F, rank.max = 100)


bm <- bm.corrected %>% 
  dplyr::rename(cn = y) %>% 
  dplyr::mutate(response = "corrected") %>% 
  dplyr::bind_rows(bm.uncorrected %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "uncorrected")) %>% 
  dplyr::bind_rows(bm.adaptive.m %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_mild")) %>% 
  dplyr::bind_rows(bm.adaptive.h %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_heavy")) %>%
  dplyr::bind_rows(bm.adaptive.mc %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_mild_complete")) %>% 
  dplyr::bind_rows(bm.adaptive.hc %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_heavy_complete"))
  
bm %<>% 
  dplyr::left_join(lfc_df %>%
                     dplyr::distinct(CompoundPlate, SampleID, cn)) %>% 
  dplyr::left_join(compound_annotations)


tr <- bm %>% 
  dplyr::distinct(SampleID, feature, GeneSymbolOfTargets) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(is.target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::group_by(SampleID, feature) %>% 
  dplyr::summarise(is.target = any(is.target, na.rm = T)) %>% 
  dplyr::ungroup()

bm %<>% 
  dplyr::left_join(tr)
  

bm %>% 
  write_csv("OncRef data/bm_lfc_comparison_all.csv")


# REST ----

file = "~/code/depmap_data/depmap_24Q4_internal.h5"
source("scripts/biomarker/biomarker_functions.R")


collapsed_l2fc <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqLFCCollapsed_adaptive_all.csv")

compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundListExtended.csv")

mut <- read_dataset(file, dataset = "Mutation")
exp <- read_dataset(file,  dataset = "Expression")
fus <- read_dataset(file,  dataset = "Fusion")
cri <- read_dataset(file,  dataset = "CRISPR")


lfc_df <- collapsed_l2fc %>%
  dplyr::inner_join(compound_annotations) 


lfc_df %<>% 
  dplyr::distinct(SampleID, CompoundPlate, CompoundName) %>% 
  dplyr::group_by(SampleID, CompoundPlate) %>% 
  dplyr::top_n(1, CompoundName) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(lfc_df)
  

air.lauc.z <- lfc_df %>%
  dplyr::distinct(CompoundName, CompoundPlate) %>%
  dplyr::rowwise() %>%
  dplyr::filter(grepl("AIR", CompoundName) | grepl("igg", CompoundName, ignore.case = T) | grepl("Fc", CompoundName, ignore.case = T)) %>%
  dplyr::filter(substr(CompoundPlate, 1,3) == "PAP") %>% 
  dplyr::left_join(lfc_df) %>%
  dplyr::mutate(fc = pmin(2^median_l2fc_uncorrected, 1)) %>%
  dplyr::filter(is.finite(fc)) %>%
  dplyr::mutate(cn = paste0(substr(CompoundName, 1, 14), "::", CompoundPlate)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "fc", fun.aggregate = mean) %>%
  log2() %>%
  scale()

air.lauc.z[is.na(air.lauc.z)] <- 0


air.lfc.z <- lfc_df %>%
  dplyr::distinct(CompoundName, CompoundPlate) %>%
  dplyr::rowwise() %>%
  dplyr::filter(grepl("AIR", CompoundName) | grepl("igg", CompoundName, ignore.case = T) | grepl("Fc", CompoundName, ignore.case = T)) %>%
  dplyr::filter(substr(CompoundPlate, 1,3) == "PAP") %>% 
  dplyr::left_join(lfc_df) %>%
  dplyr::mutate(lfc = median_l2fc) %>% # ! 
  dplyr::filter(is.finite(lfc)) %>%
  dplyr::mutate(cn = paste0(substr(CompoundName, 1, 14), "::", pert_dose, "::",  CompoundPlate)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "lfc", fun.aggregate = mean) %>%
  scale()

air.lauc.z[is.na(air.lauc.z)] <- 0
air.lfc.z[is.na(air.lfc.z)] <- 0


DELPA <- lfc_df %>%
  dplyr::distinct(CompoundName, CompoundPlate) %>%
  dplyr::rowwise() %>%
  dplyr::filter(grepl("AIR", CompoundName) | grepl("igg", CompoundName, ignore.case = T) | grepl("Fc", CompoundName, ignore.case = T)) %>%
  dplyr::filter(substr(CompoundPlate, 1,3) == "PAP") %>% 
  dplyr::left_join(lfc_df) %>%
  dplyr::filter(substr(CompoundName,1,4) == "DELP") %>% 
  dplyr::mutate(fc = pmin(2^median_l2fc_uncorrected, 1)) %>%
  dplyr::filter(is.finite(fc)) %>% 
  reshape2::acast(depmap_id ~ pert_dose, value.var = "fc", fun.aggregate = mean) 




# IGG
zz <- apply(air.lauc.z[, grepl("IGG", colnames(air.lauc.z))] < -2, 1, sum)
# zz <- apply(air.lauc.z, 1, median)

# zz3 <- air.lauc[, word(colnames(air.lauc.z)) == "DELPACIBART"]
zz3 <- rowMeans(DELPA > 0.5, na.rm = T) 
zz2.r <- rowMeans(air.lauc.z > 0)
zz2.s <- rowMeans(air.lauc.z < -2 )


cl = intersect(rownames(air.lauc.z), rownames(exp)) %>% 
  intersect(rownames(DELPA)) %>% 
  intersect(rownames(mut))



cc = ifelse(grepl("IGG", colnames(air.lauc.z)) | grepl("IgG1k", colnames(air.lauc.z)), "NC", 
            ifelse(grepl("DELP", colnames(air.lauc.z)), "PC", "Other"))


m <- scale(exp[cl, c("CCDC88C", "CD63")])
cl <- names(sort(m[,1] - m[,2]))

gr <- filtered_counts %>% 
  dplyr::distinct(depmap_id, growth_condition) %>% 
  dplyr::filter(depmap_id %in% cl) %>% 
  dplyr::mutate(growth_condition = ifelse(growth_condition == "S&A", "Mixed", growth_condition)) %>% 
  dplyr::group_by(depmap_id) %>%
  dplyr::summarise(gr = paste0(sort(unique(growth_condition)), collapse = ",")) %>% 
  dplyr::ungroup()
  


g <- gr$gr
names(g) <- gr$depmap_id


scores <- data.table::fread("~/Downloads/scores")


col_anno <- ComplexHeatmap::HeatmapAnnotation(
  #OR4F6.Fus = fus[cl, "OR4F6--"], 
  GR = g[cl],
  #DELPA.AUC =  rowMeans(DELPA[cl,], na.rm = T),
  #Sens.Frac = zz2.s[cl], 
  #Res.Frac = zz2.r[cl],
  #SLC46A3 = exp[cl, "SLC46A3"],
  CD63 = scale(exp[cl, "CD63"]),
  CCDC88C = scale(exp[cl, "CCDC88C"]),
  GNA13.DAM = mut[cl, "GNA13.DAM"])



H1 <- ComplexHeatmap::Heatmap(pmax(pmin(t(air.lauc.z[cl, ]),3),-3), 
                              show_column_names = F,
                              show_row_names = T,
                        row_split = cc,
                         column_split = ifelse(zz[cl] > 1,
                                               "All\nSens.",
                                               ifelse( zz3[cl] > .5  &   zz2.r[cl] > 0.8
                                                       , "All\nRes.", "Others")),
                        top_annotation = col_anno,
                        name = "Z-scored\nlog2(AUC)",
                        cluster_columns = F
                        ) 


H1

H2 <- ComplexHeatmap::Heatmap(pmax(pmin(t(air.lfc.z[cl, ]),3),-3), 
                              show_column_names = F, show_row_names = F,
                              column_split = ifelse(zz[cl] > 1,
                                                    "All\nSens.",
                                                    ifelse( zz3[cl] > .5  &   zz2.r[cl] > 0.8
                                                            , "All\nRes.", "Others")),
                              top_annotation = col_anno,
                        name = "Z-scored\nlog2(FC)",
                        cluster_columns = F
) 
                          
H2


ht_list = H2 %v% H1
draw(ht_list)

red <- tibble(status = ifelse(zz[cl] > 1,
                              "All\nSens.",
                              ifelse( zz3[cl] > .5  &   zz2.r[cl] > 0.8
                                      , "All\nRes.", "Others")),
       depmap_id = cl) 


tibble(depmap_id = cl, mut = mut[cl, "GNA13.DAM"]) %>% 
  dplyr::left_join(red) %>% 
  dplyr::count(mut, status)


red %>% 
  dplyr::count(status)

lfc_df %<>%
  dplyr::distinct(CompoundName, CompoundPlate) %>%
  dplyr::rowwise() %>%
  dplyr::filter(grepl("AIR", CompoundName) | grepl("igg", CompoundName, ignore.case = T) | grepl("Fc", CompoundName, ignore.case = T)) %>%
  dplyr::filter(substr(CompoundPlate, 1,3) == "PAP") %>% 
  dplyr::left_join(lfc_df) %>%
  dplyr::left_join(red) 
  
lfc_df %<>% 
  dplyr::mutate(cn = paste0(CompoundPlate, "::", SampleID, "::", pert_dose))


bm.corrected.air.full <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.corrected.air.no.S <- lfc_df %>% 
  dplyr::filter(status != "All\nSens." | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)





bm.corrected.air.no.SR <- lfc_df %>% 
  dplyr::filter(status == "Others" | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.uncorrected.air.full <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_uncorrected)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_uncorrected", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.uncorrected.air.no.S <- lfc_df %>% 
  dplyr::filter(status != "All\nSens." | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_uncorrected)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_uncorrected", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.uncorrected.air.no.SR <- lfc_df %>% 
  dplyr::filter(status == "Others" | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_uncorrected)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_uncorrected", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file,v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)




bm.adaptive.m.air.full <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_mild)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_mild", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.m.air.no.S <- lfc_df %>% 
  dplyr::filter(status != "All\nSens." | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_mild)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_mild", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.m.air.no.SR <- lfc_df %>% 
  dplyr::filter(status == "Others" | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_mild)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_mild", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)




bm.adaptive.h.air.full <- lfc_df %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_heavy)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_heavy", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.h.air.no.S <- lfc_df %>% 
  dplyr::filter(status != "All\nSens." | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_heavy)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_heavy", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100)


bm.adaptive.h.air.no.SR <- lfc_df %>% 
  dplyr::filter(status == "Others" | is.na(status)) %>% 
  dplyr::filter(is.finite(median_l2fc_adaptive_heavy)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "median_l2fc_adaptive_heavy", fun.aggregate = median) %>% 
  univariate_biomarker_table(file = file, v.X.min = 0.005,
                             features = c("CRISPR" , "CopyNumber" , "Expression", "Fusion", "Mutation", "RNAi"),
                             regression_coef = F, stability_score = F, rank.max = 100, )








bm <- bm.corrected.air.full %>% 
  dplyr::rename(cn = y) %>% 
  dplyr::mutate(response = "corrected", type = "full") %>% 
  dplyr::bind_rows(bm.corrected.air.no.S %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "corrected", type = "no.S")) %>% 
  dplyr::bind_rows(bm.corrected.air.no.SR %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "corrected", type = "no.SR")) %>% 
  dplyr::bind_rows(bm.uncorrected.air.full %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "uncorrected", type = "full")) %>% 
  dplyr::bind_rows(bm.uncorrected.air.no.S %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "uncorrected", type = "no.S")) %>% 
  dplyr::bind_rows(bm.uncorrected.air.no.SR %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "uncorrected", type = "no.SR"))  %>% 
  dplyr::bind_rows(bm.adaptive.m.air.full %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_mild", type = "full")) %>% 
  dplyr::bind_rows(bm.adaptive.m.air.no.S %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_mild", type = "no.S")) %>% 
  dplyr::bind_rows(bm.adaptive.m.air.no.SR %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_mild", type = "no.SR"))  %>% 
  dplyr::bind_rows(bm.adaptive.h.air.full %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_heavy", type = "full")) %>% 
  dplyr::bind_rows(bm.adaptive.h.air.no.S %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_heavy", type = "no.S")) %>% 
  dplyr::bind_rows(bm.adaptive.h.air.no.SR %>% 
                     dplyr::rename(cn = y) %>% 
                     dplyr::mutate(response = "adaptive_heavy", type = "no.SR")) 
  


bm %<>% 
  dplyr::left_join(lfc_df %>%
                     dplyr::distinct(CompoundPlate, SampleID, cn)) %>% 
  dplyr::left_join(compound_annotations)


tr <- bm %>% 
  dplyr::distinct(SampleID, feature, GeneSymbolOfTargets) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(is.target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::group_by(SampleID, feature) %>% 
  dplyr::summarise(is.target = any(is.target, na.rm = T)) %>% 
  dplyr::ungroup()

bm %<>% 
  dplyr::left_join(tr)


bm %>% 
  write_csv("OncRef data/bm_lfc_comparison_air_filtered.csv")

red %>% 
  write_csv("OncRef data/air/cell_line_status.csv")


red %>% 
  dplyr::count(status)
