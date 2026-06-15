library(tidyverse)
library(magrittr)
source("scripts/utils/kitchen_utensils.R")

# -----
# Load the raw data ----
# ----

compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundList26Q3.csv") %>% 
  dplyr::filter(Readout == "Sequencing-AIR")


filtered_counts <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqFilteredCounts26Q3.csv") %>%
  dplyr::filter(CompoundPlate %in% compound_annotations$CompoundPlate)


# ----
# Normalize relative to the spike-in control barcodes ----
# ----

source("scripts/normalize/normalize_functions.R")

# Be careful ! 
cb_meta <- filtered_counts %>%
  dplyr::distinct(cb_ladder, cb_name, cb_log2_dose) %>%
  tidyr::drop_na() %>% 
  dplyr::group_by(cb_ladder, cb_name) %>% 
  dplyr::summarize(cb_log2_dose = median(cb_log2_dose)) %>% 
  dplyr::ungroup()


# temp fix
cb_annots <- flag_control_bcs(filtered_counts, cb_meta, 
                              id_cols = c("pcr_plate", "pcr_well"), negcon_type = "ctl_vehicle",
                              cb_mad_cutoff = 1, req_negcon_reps = 16) %>% 
  dplyr::group_by(pcr_plate, pcr_well) %>% 
  dplyr::mutate(keep_cb = ifelse(median(n, na.rm = T) > 40, keep_cb, "No"))


normalized_counts = normalize(filtered_counts, 
                              cb_annots = cb_annots, 
                              id_cols = c("pcr_plate", "pcr_well"),
                              pseudocount = 0)


normalized_counts = add_pseudovalue(normalized_counts, 
                                    negcon_cols = c("pcr_plate", "pert_vehicle"),
                                    read_detection_limit = 10, negcon_type = "ctl_vehicle")

# -----
# QC 
# -----

source("scripts/well_metrics/well_metrics_functions.R")


cb_metrics <- calculate_cb_metrics(normalized_counts, cb_meta,
                                   group_cols = c("pcr_plate", "pcr_well"),
                                   pseudocount = 0)

normalized_counts <- normalized_counts %>% 
  dplyr::semi_join(cb_metrics %>% 
                     dplyr::filter(cb_mae < 1, cb_spearman > 0.8),
                   by = c("pcr_plate", "pcr_well"))


  

# FILTER NEGATIVE CONTROL WELLS THAT CONTROL BARCODES DOESN'T COVER AT LEAST 25% OF THE RANGE OF CELL LINE BARCODES OR
# AND ANY WELL THAT HAS SPEARMAN CORRELATION LESS THAN 0.85 WITH INTENDED DOSES. 
# CB.quantiles <- normalized_counts %>% 
#   dplyr::group_by(pcr_plate, pcr_well) %>% 
#   dplyr::arrange(n) %>% 
#   dplyr::mutate(quant = (1:n()) / n()) %>% 
#   dplyr::filter(!is.na(cb_name)) %>% 
#   dplyr::distinct(pcr_plate, pcr_well, pert_type, lua, n, quant) %>%
#   dplyr::group_by(pcr_plate, pcr_well) %>% 
#   dplyr::mutate(range = max(quant, na.rm = T) - min(quant, na.rm = T), n.b = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter((range >= 0.25) | (pert_type != "ctl_vehicle"), n.b > 5) 
# 
# 
# # # Which wells are removed
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
# 
# 
# normalized_counts <- CB.quantiles %>% 
#   dplyr::left_join(normalized_counts) %>% 
#   dplyr::group_by(pcr_plate, pcr_well, pert_type) %>%
#   dplyr::summarise(concordance = cor(n, cb_log2_dose, use = "p", method = "spearman")) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::filter(concordance >= 0.85) %>% 
#   dplyr::distinct(pcr_plate, pcr_well) %>%
#   dplyr::left_join(normalized_counts)
# 
# rm(CB.quantiles) 


source("scripts/qc_tables/qc_tables_functions.R")


normalized_counts_rm_cbc = filter_control_barcodes(normalized_counts)

outlier_pools <- normalized_counts %>% 
  filter_control_barcodes %>% 
  get_outlier_pools(ctrl_types = c("ctl_vehicle", "trt_poscon"),
                    id_cols = c("pcr_plate", "pcr_well"),
                    pcr_plate_col = "pcr_plate",
                    cell_line_cols = c("cell_set", "pool_id", "depmap_id", "lua"),
                    pool_cols = c("cell_set", "pool_id"),
                    negcon_cols = c("pcr_plate", "pert_vehicle")) %>% 
  dplyr::filter(!is.na(plate_flag))


normalized_counts <- normalized_counts %>% 
  filter_control_barcodes() %>% 
  dplyr::anti_join(outlier_pools, by = c("pcr_plate", "pcr_well", "pool_id"))


# Plate cell line QCs ----


# Generate plate cell QCs
qc_table = generate_plate_cell_table(
  normalized_counts = normalized_counts,
  ctrl_cell_line_cols = c("CompoundPlate", "cell_set", "pool_id", "depmap_id", "lua","pcr_plate", "pert_vehicle" ),  
  cell_line_cols = c("cell_set", "pool_id", "depmap_id", "lua"),
  sig_cols = c("SampleID", "pert_dose"),
  pert_plate_col = "CompoundPlate",
  pseudocount = 0,
  contains_poscon = TRUE,
  poscon = "trt_poscon", negcon = "ctl_vehicle",
  nc_variability_threshold = 1,
  error_rate_threshold = 0.05,
  pc_viability_threshold = 0.25,
  nc_raw_count_threshold = 40
)


# Add flags to plate cell QC table
qc_table = plate_cell_qc_flags(
  plate_cell_table = qc_table,
  nc_variability_threshold = 1,
  error_rate_threshold = 0.05,
  pc_viability_threshold = 0.25,
  nc_raw_count_threshold = 40,
  contains_poscon = T,
  poscon = "trt_poscon", negcon = "ctl_vehicle"
)

# TODO : choose relevant columns 
# For AIR plates: filter out suspicious cell lines 



# ----
# Compute log2-fold-changes & correct for growth condition and negative control abundance.
# ----

source("scripts/compute_l2fc/compute_l2fc_functions.R")
source("scripts/bias_correction/bias_correction_functions.R")
source("scripts/collapse_replicates/collapse_replicates_functions.R")


l2fc = compute_l2fc(normalized_counts = normalized_counts,
                    control_type = "ctl_vehicle",
                    sig_cols = c("CompoundPlate", "SampleID","pert_dose", "pert_dose_unit"),
                    ctrl_cols = c("cell_set", "pcr_plate", "replicate_plate"),
                    log2_norm_col = "log2_normalized_n",
                    cell_line_col = c("cell_set", "lua", "depmap_id", "pool_id", "growth_condition")) %>% 
  dplyr::semi_join(qc_table %>% 
                     dplyr::filter(qc_pass_pert_plate)) %>% 
  dplyr::mutate(negcon_log2_norm_n = log2(control_median_normalized_n))


l2fc <- l2fc %>% 
  dplyr::group_split(bio_rep, CompoundPlate, SampleID, pert_dose, pert_dose_unit) %>%
  lapply(apply_bias_correction, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_condition") %>%
  dplyr::bind_rows() %>% 
  dplyr::ungroup()




# Collapse l2fcs ----

collapsed_l2fc = collapse_bio_reps(l2fc = l2fc, 
                                   sig_cols = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit"), 
                                   cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"))


# -----
# Identify outlier pools and remove from l2fc tables  -----
# -----


failed_pools = get_monotonicity(collapsed_l2fc, 
                                trt_cl_cols = c("CompoundPlate", "SampleID", "depmap_id", "cell_set", "pool_id", "lua")) %>%
  flag_breaks(trt_pool_cols = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "cell_set", "pool_id"),
              trt_cell_set_cols = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "cell_set"), 
              pool_cutoff = 0.25) %>%
  dplyr::filter(!is.na(cell_set_flag))

collapsed_l2fc = collapsed_l2fc %>%
  dplyr::anti_join(failed_pools, 
                   by = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit", "cell_set", "pool_id"))


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
  write_csv("OncRef data/drc_table_AIR.csv")



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
  write_csv("OncRef data/PRISMOncologyReferenceSeqQC_AIR.csv")

l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFC_AIR.csv")

collapsed_l2fc %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqLFCCollapsed_AIR.csv")

drc_table %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqDRC_AIR.csv")


# IDENTIFY BAD CELL LINES ----


drc_table %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::inner_join(compound_annotations %>% 
                      dplyr::distinct(CompoundName, CompoundPlate, SampleID)) %>% 
  dplyr::mutate(air_cont = ifelse(CompoundName ==  "DELPACIBART AIR-NHS-MMAF", "PosCon",
                                  ifelse(grepl("IGG", CompoundName), "NegCon", "Other"))) %>% 
  dplyr::group_by(depmap_id, cell_set, air_cont) %>% 
  dplyr::summarise(m = median(auc, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = air_cont, values_from = m) %>% 
  dplyr::mutate(pass = (PosCon < 6/8) & NegCon > 7/8) %>% 
  ggplot()+
  geom_point(aes(x = NegCon, y = PosCon, 
                 color = pass)) +
  geom_abline() +
  geom_abline(intercept = -1, lty = 2) 


celllines_to_keep <- drc_table %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::inner_join(compound_annotations %>% 
                      dplyr::distinct(CompoundName, CompoundPlate, SampleID)) %>% 
  dplyr::mutate(air_cont = ifelse(CompoundName ==  "DELPACIBART AIR-NHS-MMAF", "PosCon",
                                  ifelse(grepl("IGG", CompoundName), "NegCon", "Other"))) %>% 
  dplyr::group_by(depmap_id, cell_set, air_cont) %>% 
  dplyr::summarise(m = median(auc, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = air_cont, values_from = m) %>% 
  dplyr::mutate(pass = (PosCon < 6/8) & NegCon > 7/8) %>% 
  dplyr::filter(pass) %>%
  dplyr::distinct(depmap_id, cell_set) 
 

   
air_annots <- filtered_counts %>% 
  dplyr::distinct(depmap_id, cell_set, CompoundPlate, SampleID) %>% 
  dplyr::anti_join(celllines_to_keep) %>%
  dplyr::mutate(note = "Remove (AIR)") 



filtered_counts <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqFilteredCounts26Q3.csv") %>%
  dplyr::left_join(air_annots)


filtered_counts %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqFilteredCounts26Q3.csv")

