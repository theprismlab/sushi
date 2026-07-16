library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggthemes)
library(mgcv)
source("scripts/utils/kitchen_utensils.R")
source("scripts/collapse_replicates/collapse_replicates_functions.R")


# ----
# LOAD THE TRAINING DATA ----
# ----

CompoundList <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceCompoundList.csv") %>%
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(length(unique(Readout)) > 1) %>%
  dplyr::ungroup()


PRISMOncologyReferenceLumLFCCollapsed <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumLFCCollapsed.csv") %>%  
  dplyr::filter(!outlier, priority == 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose,
                  depmap_id, pool_id, cellset, LFC) %>%  
  dplyr::inner_join(CompoundList %>% 
                      dplyr::distinct(CompoundPlate, SampleID, CompoundName))

PRISMOncologyReferenceSeqLFCCollapsed <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqLFCCollapsed.csv") %>% 
  dplyr::rename(LFC = median_l2fc, cellset = cell_set) %>% 
  dplyr::filter(num_bio_reps > 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose,
                  depmap_id, pool_id, cellset, LFC) %>%  
  dplyr::inner_join(CompoundList %>% 
                      dplyr::distinct(CompoundPlate, SampleID, CompoundName))

PRISMOncologyReferenceLumQCTable <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumQCTable.csv") %>% 
  dplyr::filter(PASS) %>% 
  dplyr::distinct(CompoundPlate, cellset, depmap_id, pool_id, NC.median) %>% 
  dplyr::group_by(CompoundPlate, cellset, depmap_id, pool_id) %>% 
  dplyr::summarise(NC.median = median(NC.median, na.rm = T)) %>% 
  dplyr::ungroup() 


# Combined collapsed LFC table

# x : Lum; y : Seq
LFC_table <- PRISMOncologyReferenceLumLFCCollapsed  %>%
  dplyr::inner_join(PRISMOncologyReferenceSeqLFCCollapsed, 
                    by = c("CompoundName", "depmap_id")) %>% 
  dplyr::filter(abs(log2(pert_dose.x / pert_dose.y)) < 0.25) %>% 
  dplyr::mutate(pert_dose = pert_dose.x / 2 + pert_dose.y / 2)


LFC_table %<>% 
  dplyr::distinct(CompoundName, CompoundPlate.x, CompoundPlate.y, pert_dose) %>%
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y) %>% 
  dplyr::arrange(pert_dose, .by_group = T) %>% 
  dplyr::mutate(dose_ix = 1:n()) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(LFC_table)



LFC_table %<>% 
  dplyr::left_join(PRISMOncologyReferenceLumQCTable,
                   by = c("cellset.x" = "cellset",
                          "pool_id.x" = "pool_id",
                          "CompoundPlate.x" = "CompoundPlate",
                          "depmap_id")) %>% 
  dplyr::group_by(cellset.x, CompoundPlate.x, pert_dose, dose_ix, CompoundName) %>% 
  dplyr::mutate(biomass = log2(1 + sum(2^NC.median * pmin(2^LFC.x,1)) / sum(2^NC.median))) %>% 
  dplyr::ungroup()



selected.compounds <- LFC_table %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y) %>% 
  dplyr::mutate(nd = length(unique(pert_dose))) %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y, nd, depmap_id, cellset.x, pool_id.x) %>%
  dplyr::summarise(auc.x = mean(pmin(2^LFC.x,1), na.rm = T), auc.y = mean(pmin(2^LFC.y,1), na.rm = T)) %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y, nd) %>%
  dplyr::summarise(r = cor(log2(auc.x), log2(auc.y), use = "p"),
                   mB = abs(mean(auc.x - auc.y))) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::top_n(1, nd) %>% 
  dplyr::top_n(1, r) %>% 
  dplyr::top_n(1, -mB) %>%
  dplyr::ungroup()

selected.compounds



# ----
# Characaterize the discrepancy -----
# ----

LFC_table %>% 
  dplyr::left_join(selected.compounds) %>%  
  dplyr::mutate(lum.bin = cut(pmin(2^LFC.x,1), seq(0, 1, 0.25))) %>% 
  dplyr::mutate(seq.bin = cut(pmin(2^LFC.y,1), seq(0, 1, 0.25) )) %>% 
  dplyr::count(lum.bin, seq.bin, dose_ix) %>% 
  dplyr::group_by(dose_ix) %>% 
  dplyr::mutate(f = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = lum.bin, y = seq.bin)) +
  geom_tile(aes( fill = log2(1+f*10000)), show.legend = F) +
  geom_text(aes(label = round(f,3))) + 
  scale_fill_viridis_c(limits = c(0, NA)) +
  theme_minimal() + 
  facet_wrap(paste0("Dose Level ", dose_ix) ~ ., ncol = 4) + 
  labs(x = "Luminex FC (binned)",  y = "Sequencing FC (binned)") 



th = 1/4

temp <- LFC_table %>% 
  dplyr::filter(is.finite(LFC.x), is.finite(LFC.y)) %>% 
  dplyr::inner_join(selected.compounds) %>% 
  dplyr::group_by(CompoundPlate.x, CompoundPlate.y, CompoundName, pert_dose, dose_ix) %>% 
  dplyr::summarise(mFC.x = mean(pmin(2^LFC.x, 1)),
                   mFC.y= mean(pmin(2^LFC.y,1)),
                   nFC.x = sum(LFC.x < -2), 
                   nFC.y = sum(LFC.y < -2)) %>% 
  dplyr::ungroup() 

temp %>% 
  ggplot(aes(x = mFC.x, y = mFC.y,
             color = dose_ix)) +
  geom_point() + 
  geom_text_repel(aes(label = ifelse(abs(mFC.x - mFC.y) >= th, CompoundName, NA)), size = 3, max.overlaps = 50) + 
  geom_abline(color = "red") +
  geom_abline(intercept = th, slope = 1, color= "red", lty = 2) +
  geom_abline(intercept = -th, slope = 1, color= "red", lty = 2) +
  theme_minimal() + scale_color_viridis_c() +
  scale_x_sqrt() + scale_y_sqrt() + 
  labs(x = "Average FC (Luminex)", y = "Average FC (Sequencing)", color = "Dose\nLevel") +
  facet_wrap(dose_ix ~ ., nrow = 2) +theme(legend.position = "none")


temp %>% 
  dplyr::mutate(flag = abs(mFC.x - mFC.y) > th) %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y) %>% 
  dplyr::summarise(mFlag = sum(flag)) %>% 
  dplyr::arrange(desc(mFlag)) %>% 
  dplyr::filter(mFlag > 2 ) %>% 
  dplyr::left_join(LFC_table) %>% 
  dplyr::left_join(temp) %>%   dplyr::mutate(flag = abs(mFC.x - mFC.y) > th, s = sign(mFC.x - mFC.y)) %>% 
  dplyr::left_join(CompoundList %>% dplyr::distinct(CompoundName, GeneSymbolOfTargets)) %>% 
  #dplyr::distinct(CompoundName)
  ggplot() +
  geom_point(aes(x = pmin(2^LFC.x, 1), y = pmin(2^LFC.y,1), color = as.factor(flag*s)), size = 0.5, show.legend = F) +
  geom_abline(color= "red") +
  # scale_x_sqrt() + scale_y_sqrt() + 
  facet_grid(dose_ix ~ GeneSymbolOfTargets + CompoundName ) + # CompoundPlate.x + CompoundPlate.y ) +
  labs(x = "Luminex FC", y = "Sequencing FC") +
  theme_minimal() + scale_color_manual(values = c("darkgreen","black", "red"))

temp %>% 
  dplyr::mutate(flag = abs(mFC.x - mFC.y) > th) %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y) %>% 
  dplyr::summarise(mFlag = sum(flag)) %>% 
  dplyr::arrange(desc(mFlag)) %>% 
  dplyr::filter(mFlag > 2) 


# instances to be filtered 
inst.redact <- temp %>% 
  dplyr::mutate(flag = abs(mFC.x - mFC.y) > th) %>% 
  dplyr::group_by(CompoundName, CompoundPlate.x, CompoundPlate.y) %>% 
  dplyr::mutate(flag = ifelse(sum(flag) > 2, TRUE, FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(flag)



LFC_table_filtered <- LFC_table %>% 
  dplyr::inner_join(selected.compounds) %>% 
  dplyr::anti_join(inst.redact) 


# next let's look at the cell lines ----

# looking at the top 4 doses 
temp <- LFC_table_filtered %>% 
  dplyr::filter(is.finite(LFC.x), is.finite(LFC.y)# , dose_ix > 5
                ) %>% 
  dplyr::group_by(depmap_id, cellset.x, cellset.y, pool_id.y) %>% 
  dplyr::summarise(mFC.x = mean(pmin(2^LFC.x,1)),
                   mFC.y= mean(pmin(2^LFC.y, 1)),
                   n = n()) %>% 
  dplyr::ungroup() 


# not as clear, for now let's keep all of them.
th = 1/4
temp %>% 
  dplyr::group_by(pool_id.y) %>% 
  dplyr::mutate(lab = ifelse(sum(abs(mFC.x - mFC.y) > th) > 1 , pool_id.y, " OTHER")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 100) %>% 
  ggplot(aes(x = mFC.x, y = mFC.y)) +
  geom_point() + 
  geom_text_repel(aes(label = ifelse(abs(mFC.x - mFC.y) > th, depmap_id, NA)),
                  show.legend = F) + 
  geom_abline(color = "red") +
  geom_abline(intercept = th, slope = 1, color= "red", lty = 2) +
  geom_abline(intercept = -th, slope = 1, color= "red", lty = 2) +
  theme_minimal() + 
  facet_wrap(cellset.x ~ .) +
  labs(x = "Average FC (Luminex)", y = "Average FC (Sequencing)", color = "Dose\nLevel") 



# ----


df <- LFC_table_filtered %>% 
  dplyr::mutate(
    LFC.x.c = pmax(pmin(LFC.x, 0), -3), 
    LFC.y.c = pmax(pmin(LFC.y, 0), -3), 
    divergence = LFC.y.c - LFC.x.c) %>% 
  tidyr::drop_na(biomass, divergence) %>%
  as.data.frame() 


seed = 25; nc = 25
set.seed(seed)
df %>%
  dplyr::distinct(depmap_id, cellset.x) %>%
  dplyr::group_by(cellset.x) %>% 
  dplyr::sample_n(nc) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(df) %>% 
  ggplot() +
  geom_point(aes(y =  pmin(2^LFC.x.c,1), # LFC.x + NC.median,
                 x = biomass, 
                 #alpha = pmin(abs(bias_raw2),3), 
                 color = divergence),
             size = 2) + #pmax(pmin(bias_raw2,3),-3))) +
  scale_color_gradient2(limits = c(-3,3)) +
  theme_minimal() + theme(legend.position = "bottom") + 
  #facet_wrap(cellset.x ~ depmap_id, ncol = 3) +
  labs(color = "LFC.seq - LFC.lum\n(capped LFC's)",
       alpha = "abs(bias_raw)\n(LFC.x and LFC.y\nare capped -3 and 0)",
       # title = "Raw", 
       x = "Total Biomass", y = "Luminex FC (capped)") 


seed = 25; nc = 5
set.seed(seed)
df %>%
  dplyr::distinct(depmap_id, cellset.x) %>%
  dplyr::group_by(cellset.x) %>% 
  dplyr::sample_n(nc) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(df) %>% 
  ggplot() +
  geom_point(aes(y =  pmin(2^LFC.x.c,1), # LFC.x + NC.median,
                 x = biomass, 
                 #alpha = pmin(abs(bias_raw2),3), 
                 color = divergence),
             size = 2) + #pmax(pmin(bias_raw2,3),-3))) +
  scale_color_gradient2(limits = c(-3,3)) +
  theme_minimal() + theme(legend.position = "bottom") + 
  facet_wrap(cellset.x ~ depmap_id, ncol = 5) +
  labs(color = "LFC.seq - LFC.lum\n(capped LFC's)",
       alpha = "abs(bias_raw)\n(LFC.x and LFC.y\nare capped -3 and 0)",
       # title = "Raw", 
       x = "Total Biomass", y = "Luminex FC (capped)") 





global_model <- gam(
  divergence ~ te(biomass, LFC.x.c,  k = c(5, 5), bs = 'ts'),
  data = df, 
  method = "REML",
  select = TRUE
)

summary(global_model)

C = 1/2


df$predicted_divergence <-  as.numeric(predict(global_model, newdata = df))
df$lfc_weight <- 1 - exp(-(df$LFC.x.c / C)^2)
df$pinched_divergence <- df$predicted_divergence * df$lfc_weight


# Visualization of the model -----

grid <- expand.grid(
  biomass = seq(min(df$biomass), max(df$biomass), length.out = 100),
  LFC.x.c = seq(-3, 0, length.out = 100)
)

grid$Penalty <- as.numeric(predict(global_model, newdata = grid))
grid$lfc_weight <- 1 - exp(-(grid$LFC.x.c / C)^2)
grid$Pinched_Penalty <- grid$Penalty * grid$lfc_weight


p1 = ggplot(grid, aes(x = biomass, y = LFC.x.c, fill = Penalty)) +
  geom_raster(alpha = 1) +
  geom_contour(aes(z = Penalty), color = "black", alpha = 0.2, bins = 15) +
  geom_hline(aes(yintercept = 0), color = "gray", lwd = 0.5, lty = 2, inherit.aes = F) + 
  # 2. The Highlighted ZERO Contour
  geom_contour(aes(z = Penalty), breaks = 0, color = "black", linewidth = 1, linetype = "dashed") +
  # The Zero-Centered Color Scale
  scale_fill_gradient2(
    low = "#B2182B",       
    mid = "white",         
    high = "#2166AC",      
    midpoint = 0,   
    limits = c(-4, 1.5), 
    name = "Effective\nPenalty"
  ) +
  
  theme_base() +
  labs(
    # title = "Correction Curves",
    subtitle = "Predicted divergence", # after shifting for negative control",
    x = "Total Biomass",
    y = "Luminex LFC (capped)"
  ) 



p2 = ggplot(grid, aes(x = biomass, y = LFC.x.c, fill = Pinched_Penalty)) +
  geom_raster(alpha = 1) +
  geom_contour(aes(z = Pinched_Penalty), color = "black", alpha = 0.2, bins = 15) +
  geom_hline(aes(yintercept = 0), color = "gray", lwd = 0.5, lty = 2, inherit.aes = F) + 
  # 2. The Highlighted ZERO Contour
  geom_contour(aes(z = Pinched_Penalty), breaks = 0, color = "black", linewidth = 1, linetype = "dashed") +
  # The Zero-Centered Color Scale
  scale_fill_gradient2(
    low = "#B2182B",       
    mid = "white",         
    high = "#2166AC",      
    midpoint = 0,   
    limits = c(-4, 1.5), 
    name = "Effective\nPenalty"
  ) +
  theme_base() +
  labs(
    # title = "Correction Surves",
    subtitle = "Predicted divergence after protecting LFC ~ 0",
    x = "Total Biomass",
    y = "Luminex LFC (capped)"
  ) 


cowplot::plot_grid(p1,p2, nrow = 1)

# ----- 
# TO-DO: Benchmarks! 
# -----


# ----
# TRANSFORMING AND PROCESSING THE REST OF THE DATA  ----
# ----

PRISMOncologyReferenceLumQCTable <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumQCTable.csv")
PRISMOncologyReferenceLumAnalyteMeta <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumAnalyteMeta.csv")
PRISMOncologyReferenceLumInstMeta <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumInstMeta.csv")


data.to.keep <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumLFCCollapsed.csv") %>%  
  dplyr::filter(!outlier, priority == 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pool_id, depmap_id, screen) %>%
  dplyr::left_join(PRISMOncologyReferenceLumAnalyteMeta %>% 
                     dplyr::distinct(depmap_id, pool_id, cellset, screen, analyte_id, note)) %>%
  dplyr::filter(is.na(note)) 



PRISMOncologyReferenceLumLFC <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceLumLFC.csv") %>% 
  dplyr::filter(is.finite(LFC)) %>% 
  dplyr::inner_join(PRISMOncologyReferenceLumInstMeta) %>% 
  dplyr::inner_join(PRISMOncologyReferenceLumAnalyteMeta) %>% 
  dplyr::inner_join(PRISMOncologyReferenceLumQCTable %>% 
                     dplyr::distinct(prism_replicate, cellset, analyte_id, pool_id, NC.median)) %>% 
  dplyr::group_by(cellset, prism_replicate, pert_well) %>% 
  dplyr::mutate(biomass = log2(1 + sum(2^NC.median * pmin(2^LFC,1)) / sum(2^NC.median)),
                LFC.x.c = pmax(pmin(LFC, 0), -3)) %>% 
  dplyr::ungroup() 


PRISMOncologyReferenceLumLFC$predicted_divergence <-  as.numeric(predict(global_model, newdata = PRISMOncologyReferenceLumLFC))
PRISMOncologyReferenceLumLFC$lfc_weight <- 1 - exp(-(PRISMOncologyReferenceLumLFC$LFC.x.c / C)^2)
PRISMOncologyReferenceLumLFC$pinched_divergence <- PRISMOncologyReferenceLumLFC$predicted_divergence * PRISMOncologyReferenceLumLFC$lfc_weight


# !!!

l2fc_lum_harmonized <- PRISMOncologyReferenceLumLFC %>% 
  dplyr::inner_join(data.to.keep) %>% 
  dplyr::mutate(l2fc_transformed = LFC  + pinched_divergence, pert_type = "trt_cp") %>% 
  dplyr::rename(pcr_plate = prism_replicate, cell_set = cellset, lua = analyte_id, replicate_plate =replicate,
                l2fc = LFC,
                l2fc_uncorrected = LFC_uncorrected) %>% 
  dplyr::distinct(CompoundPlate, pcr_plate, SampleID, pert_dose, pert_dose_unit, 
                  cell_set, pool_id, lua, depmap_id, pert_type, replicate_plate, 
                  l2fc_transformed, l2fc, l2fc_uncorrected)


collapsed_l2fc_lum_harmonized = collapse_bio_reps(l2fc = l2fc_lum_harmonized, 
                                                   median_cols = c("l2fc", "l2fc_uncorrected", "l2fc_transformed"),
                                   sig_cols = c("CompoundPlate", "SampleID", "pert_dose", "pert_dose_unit"), 
                                   cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"))



source("scripts/drc/dose_response_functions.R")


plates <- l2fc_lum_harmonized %>% 
  dplyr::distinct(CompoundPlate, SampleID) %>% 
  dplyr::count(CompoundPlate) %>% 
  dplyr::arrange(n) %>%
  .$CompoundPlate

drc_table_lum_harmonized <- list(); 
ix = 1;   now = Sys.time()

for(plate in plates){
  print(plate)
  print(Sys.time() - now)
  now = Sys.time()
  
  drc_table_lum_harmonized[[ix]] <- l2fc_lum_harmonized %>% 
    dplyr::filter(CompoundPlate == plate) %>% 
    dplyr::group_split(SampleID) %>% 
    parallel::mclapply(create_drc_table,
                       l2fc_col = "l2fc_transformed",
                       screen_type = "MTS_SEQ", cell_line_cols = c("depmap_id", "lua", "pool_id", "cell_set"),
                       treatment_cols = c("CompoundPlate", "SampleID"), dose_col = "pert_dose", type_col = "pert_type", cap_for_viability = 1.5,
                       mc.cores = pmax(parallel::detectCores() - 2,1)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(successful_fit)
  
  print(ix)
  print(plate)
  print(Sys.time() - now)
  now = Sys.time()
  

  ix <- ix + 1 
}


drc_table_lum_harmonized <- drc_table_lum_harmonized %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(response = "harmonized_luminex") 




# Compute fitted l2fc values ----

fitted_l2fc <- collapsed_l2fc_lum_harmonized %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, pert_dose_unit) %>% 
  dplyr::inner_join(drc_table_lum_harmonized) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(l2fc_fitted = log2(lower_limit + (upper_limit - lower_limit) / (1 + (pert_dose / inflection)^-slope))) %>%  
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, pert_dose_unit, cell_set, depmap_id, pool_id, lua, l2fc_fitted, response) 


collapsed_l2fc_lum_harmonized %<>% 
  dplyr::left_join(fitted_l2fc) 



# SAVE FILES ----

l2fc_lum_harmonized %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLumLFC.csv")
  
collapsed_l2fc_lum_harmonized %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLumLFCCollapsed.csv")
  

drc_table_lum_harmonized %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLumDRC.csv")


# Portal files --- !

compound_annotations <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceCompoundList.csv")  %>% 
  dplyr::filter(Prioritized, hasData) %>% 
  dplyr::group_by(CompoundName) %>%
  dplyr::top_n(1, Readout) %>% 
  dplyr::ungroup() 
  
  

seqResponseCurves <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqResponseCurves.csv") %>% 
  dplyr::mutate(EC50 = as.numeric(EC50))

drc_table_lum_harmonized %>% 
  dplyr::inner_join(compound_annotations) %>% 
  dplyr::filter(successful_fit) %>% 
  dplyr::rename(ModelID = depmap_id, EC50 = inflection, LowerAsymptote = lower_limit, UpperAsymptote = upper_limit, Slope = slope) %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  dplyr::bind_rows(seqResponseCurves) %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedResponseCurves.csv")



seqLog2AUC <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqLog2AUCMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix()


drc_table_lum_harmonized %>% 
  dplyr::inner_join(compound_annotations) %>% 
  dplyr::filter(successful_fit) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2_auc") %>% 
  reshape2::melt() %>% 
  dplyr::bind_rows(reshape2::melt(seqLog2AUC)) %>% 
  dplyr::filter(is.finite(value)) %>% 
  reshape2::acast(Var1 ~ Var2, value.var = "value") %>% 
  write.csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLog2AUCMatrix.csv")




seqLFCCollapsed <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqLog2ViabilityCollapsedMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix()


collapsed_l2fc_lum_harmonized %>%
  dplyr::inner_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate)) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, depmap_id, l2fc_fitted) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "l2fc_fitted") %>% 
  reshape2::melt() %>% 
  dplyr::bind_rows(reshape2::melt(seqLFCCollapsed)) %>% 
  dplyr::filter(is.finite(value)) %>% 
  reshape2::acast(Var1 ~ Var2, value.var = "value") %>% 
  write.csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLog2ViabilityCollapsedMatrix.csv")



seqLFCCollapsedConditions <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqLog2ViabilityCollapsedConditions.csv")



collapsed_l2fc_lum_harmonized %>%
  dplyr::inner_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate)) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, l2fc_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit) %>% 
  dplyr::bind_rows(seqLFCCollapsedConditions) %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedLog2ViabilityCollapsedConditions.csv")




LFC_ <- compound_annotations %>% 
  dplyr::distinct(SampleID, CompoundName, CompoundPlate) %>% 
  dplyr::inner_join(l2fc_lum_harmonized) %>% 
  dplyr::filter(is.finite(l2fc_transformed)) %>% 
  dplyr::mutate(l2fc = l2fc_transformed) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate, l2fc) %>% 
  dplyr::group_by(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::mutate(Replicate = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) 


seqViabilityConditions <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqViabilityConditions.csv")

seqViabilityMatrix <- data.table::fread("OncRef data/release files/PRISMOncologyReferenceSeqViabilityMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix()

colnames(seqViabilityMatrix) <- seqViabilityMatrix[1, ] 
seqViabilityMatrix <- seqViabilityMatrix[-1,] 


seqLFC_ <- seqViabilityMatrix%>%
  reshape2::melt() %>% 
  rename(depmap_id = Var1, Label = Var2, viability = value) %>% 
  dplyr::filter(is.finite(viability)) %>% 
  dplyr::inner_join(seqViabilityConditions)  %>% 
  dplyr::mutate(l2fc = log2(viability)) %>% 
  dplyr::select(-viability, -Label)





LFC.conditions <- LFC_ %>% 
  dplyr::bind_rows(seqLFC_) %>% 
  dplyr::distinct(SampleID, Dose, DoseUnit, CompoundPlate, Replicate) %>% 
  dplyr::arrange(SampleID, CompoundPlate, Dose, DoseUnit, Replicate) %>% 
  dplyr::mutate(Label = 1:n() - 1) 



LFC_ %>% 
  dplyr::bind_rows(seqLFC_) %>% 
  dplyr::left_join(LFC.conditions) %>% 
  dplyr::mutate(viability = pmin(2^l2fc, 1.5)) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "viability") %>%
  write.csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedViabilityMatrix.csv")

LFC.conditions %>% 
  write_csv("OncRef data/release files/PRISMOncologyReferenceHarmonizedViabilityConditions.csv")





