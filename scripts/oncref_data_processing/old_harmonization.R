library(renv)
renv::restore()
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggthemes)
library(mgcv)

source("scripts/UTILITIES.R")


# ----
# LOAD THE TRAINING DATA ----
# ----


CompoundList <- data.table::fread("data/sequencing data/PRISMOncologyReferenceCompoundList.csv") %>%
  # dplyr::distinct() %>% 
  #  dplyr::mutate(cn = paste0(SampleID, "::", CompoundPlate),
  #                cn2 = paste0(CompoundName, "::", Readout)) %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(length(unique(Readout)) > 1) %>%
  dplyr::ungroup()


PRISMOncologyReferenceLumLFCCollapsed <- data.table::fread("data/processed data/PRISMOncologyReferenceLumLFCCollapsed.csv") %>%  
  dplyr::filter(!outlier, priority == 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose,
                  depmap_id, pool_id, cellset, LFC) %>%  
  dplyr::inner_join(CompoundList %>% 
                      dplyr::distinct(CompoundPlate, SampleID, CompoundName))



PRISMOncologyReferenceSeqLFCCollapsed <- data.table::fread("data/sequencing data/PRISMOncologyReferenceSeqLFCCollapsed.csv") %>% 
  dplyr::rename(LFC = median_l2fc, cellset = cell_set) %>% 
  dplyr::filter(num_bio_reps > 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose,
                  depmap_id, pool_id, cellset, LFC) %>%  
  dplyr::inner_join(CompoundList %>% 
                      dplyr::distinct(CompoundPlate, SampleID, CompoundName))

PRISMOncologyReferenceLumQCTable <- data.table::fread("data/processed data/PRISMOncologyReferenceLumQCTable.csv") %>% 
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
# 1. Characaterize the discrepancy -----
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
  dplyr::filter(is.finite(LFC.x), is.finite(LFC.y), dose_ix > 5) %>% 
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

p1


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

df %<>% 
  dplyr::mutate(LFC.x.corrected = LFC.x + pinched_divergence) %>% 
  dplyr::mutate(corrected_divergence = LFC.y.c - pmax(pmin(LFC.x.corrected, 0), -3))


df %>% 
  dplyr::mutate(
    LFC.x.bin = cut(LFC.x.c, c(-3, -2, -1, 0), include.lowest = T) ,
    biomass.bin = cut(biomass, quantile(biomass), include.lowest = T) ) %>%
  dplyr::group_by(LFC.x.bin, biomass.bin) %>%
  dplyr::sample_frac(0.1) %>% 
  ggplot(aes(x = abs(divergence), y = abs(corrected_divergence)),) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_abline() +
  facet_grid(LFC.x.bin ~ biomass.bin)

temp <- df %>% 
  dplyr::mutate(
    LFC.x.bin = cut(LFC.x.c, seq(-3, 0, 1), include.lowest = T) ,
    LFC.y.bin = cut(LFC.y.c, seq(-3, 0, 1), include.lowest = T) ,
    biomass.bin = cut(biomass, c(0, 0.5, 0.75 ,1), include.lowest = T) ) %>%
  dplyr::filter(is.finite(divergence), is.finite(corrected_divergence), is.finite(LFC.x.c)) %>%
  dplyr::group_by(LFC.x.bin,LFC.y.bin,  biomass.bin) %>% 
  dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup()


temp %>% 
  ggplot(aes(x = biomass.bin, y = abs(divergence) - abs(corrected_divergence), color = biomass.bin)) +
  geom_boxplot(size = 0.5, alpha = 0.5) +
  geom_hline(aes(yintercept = 0), color = "black") + 
  geom_text(data = dplyr::distinct(temp, nn, LFC.x.bin, LFC.y.bin, biomass.bin, nn),
            aes(label = round(nn/sum(nn) * 100, 2), y = 2), color = "black") +
  facet_grid(paste0("Lum: ", LFC.x.bin) ~ 
               paste0("Seq: ", LFC.y.bin)) +
  theme_bw()




temp <- df %>% 
  dplyr::mutate(
    LFC.x.bin = cut(LFC.x.c, seq(-3, 0, 0.75), include.lowest = T) ,
    biomass.bin = cut(biomass, c(0,  0.5, 0.75 ,1), include.lowest = T) ) %>%
  dplyr::filter(is.finite(divergence), is.finite(corrected_divergence), is.finite(LFC.x.c)) %>% 
  dplyr::group_by(LFC.x.bin, biomass.bin) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(cols = c(divergence, corrected_divergence),
                      names_to = "cor", values_to = "divergence") 


temp %>% 
  ggplot(aes(x = LFC.x.bin, y = divergence, color = cor)) +
  geom_boxplot(size = 0.5, alpha = 0.5) +
  geom_text(data = dplyr::distinct(temp, nn, LFC.x.bin, biomass.bin, nn),
            aes(label = paste0(round(nn/sum(nn) * 100, 2), "%"), y = 2), color ="black") + 
  facet_grid(. ~ reorder(paste0("BM : ", biomass.bin), as.numeric(biomass.bin)), scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(y = "Divergence (Sequencing - Luminex LFC)",
       x = "LFC Luminex (capped)",
       color =  NULL)





p1 = temp %>% 
  tidyr::pivot_wider(names_from = "cor", values_from = "divergence") %>% 
  dplyr::group_by(depmap_id) %>%
  dplyr::summarise(corrected_divergence = mean(corrected_divergence, na.rm = T),
                   divergence = mean(divergence, na.rm = T),
                   n = n()) %>% 
  ggplot() +
  geom_point(aes(x = abs(divergence), y = abs(corrected_divergence),
                 color = n), size = 1, alpha = 0.5) +
  geom_abline() +
  scale_color_viridis_b(transform = "log2") +
  theme_bw() +
  labs(x = "Absolute Divergence (Seq - Lum LFC)", y = "Absolute Divergence Post-Correction")



p2 = temp %>% 
  tidyr::pivot_wider(names_from = "cor", values_from = "divergence") %>% 
  dplyr::mutate(sens.y = ifelse(LFC.y < -1, "Seq LFC < -1", "Seq LFC > -1"), 
                sens.x = ifelse(LFC.x < -1, "Lum LFC < -1", "Lum LFC > -1")) %>% 
  dplyr::group_by(depmap_id,sens.y, sens.x) %>%
  dplyr::summarise(corrected_divergence = mean(corrected_divergence, na.rm = T),
                   divergence = mean(divergence, na.rm = T),
                   n = n()) %>% 
  ggplot() +
  geom_point(aes(x = abs(divergence), y = abs(corrected_divergence),
                 color = n), size = 1, alpha = 0.5) +
  geom_abline() +
  # geom_text(aes(x = 0, y = 0, label = n), color = "red") + 
  facet_grid(sens.x ~ reorder(sens.y, -as.numeric(sens.y))) +
  scale_color_viridis_b(transform = "log2") +
  theme_bw() +
  labs(x = "Absolute Divergence (Seq - Lum LFC)", y = "Absolute Divergence Post-Correction")


cowplot::plot_grid(p1, p2, nrow = 1)






temp <- df %>% 
  dplyr::mutate(
    LFC.y.bin = cut(LFC.y.c, seq(-3, 0, 0.5), include.lowest = T) ,
    biomass.bin = cut(biomass, c(0, 0.5, 0.75 ,1), include.lowest = T) ) %>%
  dplyr::filter(is.finite(divergence), is.finite(corrected_divergence), is.finite(LFC.y.c)) %>% 
  dplyr::group_by(LFC.y.bin, biomass.bin) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(cols = c(divergence, corrected_divergence),
                      names_to = "cor", values_to = "divergence") 


temp %>% 
  ggplot(aes(x = LFC.y.bin, y = divergence, color = cor)) +
  geom_text(data = dplyr::distinct(temp, nn, LFC.y.bin, biomass.bin, nn),
            aes(label = paste0(round(nn/sum(nn) * 100, 2), "%"), y = 2), color ="black") + 
  geom_boxplot(size = 0.5, alpha = 0.5) +
  geom_text(data = dplyr::distinct(temp, nn, LFC.y.bin, biomass.bin, nn),
            aes(label = paste0(round(nn/sum(nn) * 100, 2), "%"), y = 2), color ="black") + 
  facet_grid(. ~ biomass.bin, scales = "free") +
  coord_flip()

df_long <- df %>%
  #dplyr::filter(LFC.y.c < -1.0) %>% 
  select(divergence, corrected_divergence, LFC.x.c, biomass, LFC.y.c) %>%
  pivot_longer(cols = c(1,2), names_to = "Metric", values_to = "Error")


ggplot(df_long, aes(x = Error, fill = Metric, color = Metric)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  
  # Ensure both fill and color scales have the exact same names, values, and labels
  scale_fill_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  scale_color_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(-3, 2)) +
  labs(
    title = "Error Distribution in Responding Cells (LFC < -1.0)",
    subtitle = "Correction eliminates the massive optical background bias in dead cells",
    x = "Divergence (Luminex - Sequencing)", 
    y = "Density"
  )  +
  facet_grid(ifelse(LFC.y.c < -1, "Seq LFC < -1", "Seq LFC > -1") ~ .)

ggplot(df_long, aes(x = Error, fill = Metric, color = Metric)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  
  # Ensure both fill and color scales have the exact same names, values, and labels
  scale_fill_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  scale_color_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(-3, 2)) +
  labs(
    x = "Divergence (Sequencing - Luminex)", 
    y = "Density"
  ) +
  facet_grid(ifelse(LFC.y.c < -1, "Seq LFC < -1", "Seq LFC > -1") ~ ~ ifelse(LFC.x.c < -1, "Luminex LFC < -1", "Luminex LFC > -1"))


ggplot(df_long, aes(x = Error, fill = Metric, color = Metric)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  
  # Ensure both fill and color scales have the exact same names, values, and labels
  scale_fill_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  scale_color_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(-3, 2)) +
  labs(
    x = "Divergence (Sequencing - Luminex)", 
    y = "Density"
  ) +
  facet_grid(ifelse(LFC.y.c < -1, "Seq LFC < -1", "Seq LFC > -1") ~ ifelse(LFC.x.c < -1, ifelse(biomass < 0.75, "Biomass < 0.75 & Lum LFC < -1", "Biomass > 0.75 & Lum LFC < -1") ,
                                                                           "Luminex LFC > -1"))



ggplot(df_long, aes(x = Error, fill = Metric, color = Metric)) +
  geom_histogram(alpha = 0.5, adjust = 1.5, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  
  # Ensure both fill and color scales have the exact same names, values, and labels
  scale_fill_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  scale_color_manual(
    name = "Metric",
    values = c("divergence" = "#00BFC4", "corrected_divergence" = "#F8766D"),
    labels = c("divergence" = "Raw Divergence", "corrected_divergence" = "Corrected Divergence")
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(-3, 2)) +
  labs(
    x = "Divergence (Sequencing - Luminex)", 
    y = "Density"
  ) +
  facet_grid(ifelse(LFC.y.c < -1, "Seq LFC < -1", "Seq LFC > -1") ~ ifelse(LFC.x.c < -1, "Lum LFC < -1" , # ifelse(biomass < 0.75, "Biomass < 0.75 & Lum LFC < -1", "Biomass > 0.75 & Lum LFC < -1") ,
                                                                           "Luminex LFC > -1"),
             scales = "free_y") 





# Plot 1: Raw Concordance
p1 <- ggplot(df, aes(x = LFC.y.c, y = LFC.x.c)) +
  geom_hex(bins = 100) + # Hexbin is much faster and cleaner for 700k points than geom_point
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Pre-Correction", x = "Sequencing LFC", y = "Raw Luminex LFC") +
  theme(legend.position = "none")


# Plot 2: Corrected Concordance
p2 <- ggplot(df, aes(x = LFC.y.c, y = pmax(pmin(LFC.x.corrected, 0), -3))) +
  geom_hex(bins = 100) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Post-Correction", x = "Sequencing LFC", y = "Corrected Luminex LFC")

cowplot::plot_grid(p1,p2)






PRISMOncologyReferenceLumQCTable <- data.table::fread("data/processed data/PRISMOncologyReferenceLumQCTable.csv")
PRISMOncologyReferenceLumAnalyteMeta <- data.table::fread("data/input data/PRISMOncologyReferenceLumAnalyteMeta.csv")
PRISMOncologyReferenceLumInstMeta <- data.table::fread("data/input data/PRISMOncologyReferenceLumInstMeta.csv")

data.to.keep <- data.table::fread("data/processed data/PRISMOncologyReferenceLumLFCCollapsed.csv") %>%  
  dplyr::filter(!outlier, priority == 1, is.finite(LFC)) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pool_id, depmap_id) 


PRISMOncologyReferenceLumLFC <- data.table::fread("data/processed data/PRISMOncologyReferenceLumLFC.csv") %>% 
  dplyr::filter(PASS, is.finite(LFC_regressed)) %>% 
  dplyr::left_join(PRISMOncologyReferenceLumQCTable %>% 
                     dplyr::distinct(prism_replicate, cellset, analyte_id, pool_id, NC.median)) %>% 
  dplyr::group_by(cellset, prism_replicate, pert_well) %>% 
  dplyr::mutate(biomass = log2(1 + sum(2^NC.median * pmin(2^LFC_regressed,1)) / sum(2^NC.median)),
                LFC.x.c = pmax(pmin(LFC_regressed, 0), -3)) %>% 
  dplyr::ungroup() %>%
  dplyr::inner_join(PRISMOncologyReferenceLumInstMeta) %>% 
  dplyr::inner_join(PRISMOncologyReferenceLumAnalyteMeta) %>% 
  dplyr::inner_join(data.to.keep) 


PRISMOncologyReferenceLumLFC$predicted_divergence <-  as.numeric(predict(global_model, newdata = PRISMOncologyReferenceLumLFC))
PRISMOncologyReferenceLumLFC$lfc_weight <- 1 - exp(-(PRISMOncologyReferenceLumLFC$LFC.x.c / C)^2)
PRISMOncologyReferenceLumLFC$pinched_divergence <- PRISMOncologyReferenceLumLFC$predicted_divergence * PRISMOncologyReferenceLumLFC$lfc_weight


PRISMOncologyReferenceLumLFC.Harmonized <- PRISMOncologyReferenceLumLFC %>% 
  dplyr::mutate(l2fc = LFC_regressed  + pinched_divergence, pert_type = "trt_cp") %>% 
  dplyr::rename(pcr_plate = prism_replicate, cell_set = cellset, lua = analyte_id, replicate_plate =replicate) %>% 
  dplyr::distinct(CompoundPlate, pcr_plate, SampleID, pert_dose, pert_dose_unit, 
                  cell_set, pool_id, lua, depmap_id, pert_type, replicate_plate, 
                  l2fc)


PRISMOncologyReferenceLumLFC.Harmonized %>% 
  dplyr::distinct(CompoundPlate, SampleID) %>% 
  dplyr::count(CompoundPlate) %>% 
  dplyr::arrange(n) %>% as.data.frame()
