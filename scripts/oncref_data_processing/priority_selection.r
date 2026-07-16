library(tidyverse)
library(magrittr)
source("scripts/utils/kitchen_utensils.R")


DRC <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqDRC26Q3.csv")

CompoundList <- data.table::fread("OncRef data/luminex_data/PRISMOncologyReferenceCompoundList.csv")

mix <- CompoundList %>% 
  dplyr::filter(Readout == "Sequencing") %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::distinct(CompoundPlate, CompoundName, SampleID, GeneSymbolOfTargets)
  



LAUC <- DRC %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::inner_join(mix) %>% 
  reshape2::acast(depmap_id ~ CompoundPlate + SampleID + maximum_dose,
                  value.var = "log2_auc", fun.aggregate = median)



source("scripts/biomarker/biomarker_functions.R")

bm <- univariate_biomarker_table(LAUC, file = "~/code/depmap_data/depmap_26Q1_internal.h5",
                          features = c("CRISPR", "Expression", "Mutation", "Fusion", "RNAi"))




bm %>% 
  dplyr::mutate(CompoundPlate = word(y,1,-3, sep = fixed("_")),
                SampleID = word(y, -2, sep = fixed("_")),
                top_dose = word(y, -1, sep = fixed("_"))) %>% 
  dplyr::left_join(mix) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = fixed(";")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(target) %>% 
  ggplot() +
  geom_point(aes(y = paste0(CompoundPlate, "/", top_dose),
                 x = rank, color = feature_set),
             size = 2, alpha = 0.25) +
  facet_wrap(CompoundName ~ ., scales = "free") +
  scale_x_sqrt()



bm %>% 
  dplyr::mutate(CompoundPlate = word(y,1,-3, sep = fixed("_")),
                SampleID = word(y, -2, sep = fixed("_")),
                top_dose = word(y, -1, sep = fixed("_"))) %>% 
  dplyr::left_join(mix) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = fixed(";")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(target) %>% 
  dplyr::filter(CompoundName == "TRASTUZUMAB EMTANSINE") %>% 
  ggplot() +
  geom_point(aes(x = abs(correlation_coef), y = rank, color = CompoundPlate,
                 size = log2(1+n)), alpha = 0.5) +
  facet_wrap(feature + feature_set ~ ., scales = "free") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA))

  

selected.seq <- CompoundList %>% 
  dplyr::filter(Readout == "Sequencing", hasData) %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(CompoundPlate, CompoundName, SampleID, GeneSymbolOfTargets, hasData) %>% 
  dplyr::distinct(CompoundName) %>% 
  dplyr::mutate(CompoundPlate = c("MTS_SEQ004_S1", "MTS_SEQ004_S1", "MTS_SEQ004_S1",
                                  "PMTS093", "MTS_SEQ004_S1", "PMTS092",
                                  "PMTS092", "MTS_SEQ004_S1", "MTS_SEQ004_S1",
                                  "MTS_SEQ004_S1", "MTS_SEQ004_S2", "PMTS090",
                                  "PMTS090", "PMTS077", "PMTS077" ,
                                  "PMTS077", "PMTS077", "PMTS077",
                                  "PAPS006")) %>% 
  dplyr::mutate(Prioritized2 = TRUE, Readout = "Sequencing")






mix <- CompoundList %>% 
  dplyr::filter(Readout == "Sequencing-AIR") %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::distinct(CompoundPlate, CompoundName, SampleID, GeneSymbolOfTargets)



LAUC <- DRC %>% 
  dplyr::filter(response == "corrected") %>% 
  dplyr::inner_join(mix) %>% 
  reshape2::acast(depmap_id ~ CompoundPlate + SampleID + maximum_dose,
                  value.var = "log2_auc", fun.aggregate = median)



bm <- univariate_biomarker_table(LAUC, file = "~/code/depmap_data/depmap_26Q1_internal.h5",
                                 features = c( "Expression"))




bm %>% 
  dplyr::mutate(CompoundPlate = word(y,1,-3, sep = fixed("_")),
                SampleID = word(y, -2, sep = fixed("_")),
                top_dose = word(y, -1, sep = fixed("_"))) %>% 
  dplyr::left_join(mix) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = fixed(";")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(target) %>% 
  ggplot() +
  geom_point(aes(y = paste0(CompoundPlate, "/", top_dose),
                 x = rank, color = feature_set),
             size = 2, alpha = 0.25) +
  facet_wrap(CompoundName ~ ., scales = "free") +
  scale_x_sqrt()


bm %>% 
  dplyr::mutate(CompoundPlate = word(y,1,-3, sep = fixed("_")),
                SampleID = word(y, -2, sep = fixed("_")),
                top_dose = word(y, -1, sep = fixed("_"))) %>% 
  dplyr::left_join(mix) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = fixed(";")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(target = grepl(GeneSymbolOfTargets, feature)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(target) %>% 
  dplyr::filter(CompoundName == "HUMAN IGG1 ISOTYPE CTRL AIR-NHS-MMAF") %>% 
  ggplot() +
  geom_point(aes(x = abs(correlation_coef), y = rank, color = CompoundPlate,
                 size = log2(1+n)), alpha = 0.5) +
  facet_wrap(feature + feature_set ~ ., scales = "free") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA))



selected.air <- CompoundList %>% 
  dplyr::filter(Readout == "Sequencing-AIR", hasData) %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(CompoundPlate, CompoundName, SampleID, GeneSymbolOfTargets, hasData) %>% 
  dplyr::left_join(tibble(n = colSums(is.finite(LAUC)) ,
                          CompoundPlate = word(colnames(LAUC), sep = fixed("_")),
                          SampleID = word(colnames(LAUC), 2, sep = fixed("_")))) %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::top_n(1, n) %>% 
  dplyr::distinct(CompoundName, CompoundPlate) %>% 
  dplyr::mutate(Prioritized2 = TRUE, Readout = "Sequencing-AIR")



CL1 <- CompoundList %>% 
  dplyr::group_by(CompoundName, Readout) %>% 
  dplyr::filter(n() == 1) %>% 
  dplyr::ungroup() 


CL2 <-CompoundList %>% 
  dplyr::filter(Readout != "Luminex") %>% 
  dplyr::group_by(CompoundName, Readout) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Prioritized = ifelse(!hasData, NA, Prioritized)) %>% 
  dplyr::left_join(selected %>% 
                     dplyr::bind_rows(selected.air)) %>% 
  dplyr::mutate(Prioritized = ifelse(is.na(Prioritized2), ifelse(hasData, FALSE, NA), TRUE)) %>% 
  dplyr::select(-Prioritized2) %>% 
  dplyr::distinct()



CL3 <- CompoundList %>% 
  dplyr::filter(Readout == "Luminex") %>% 
  dplyr::group_by(CompoundName, Readout) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Prioritized = ifelse(!hasData, NA, Prioritized)) 



CL <- CL1 %>% 
  dplyr::bind_rows(CL2) %>% 
  dplyr::bind_rows(CL3)


CL %>% 
  dplyr::arrange(CompoundName, Prioritized, Readout) %>% 
  write_csv("OncRef data/PRISMOncologyReferenceCompoundList26Q3.csv")

