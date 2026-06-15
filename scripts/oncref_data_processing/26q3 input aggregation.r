library(tidyverse)
library(useful)
library(magrittr)
library(taigr)


comprehensive_list <- data.table::fread("~/Downloads/INTERNAL_Oncology Reference Screen Info for DMC Partners - Comprehensive Luminex.csv") %>% 
  dplyr::mutate(Readout = "Luminex") %>% 
  dplyr::bind_rows(data.table::fread("~/Downloads/INTERNAL_Oncology Reference Screen Info for DMC Partners - Comprehensive Sequencing (1).csv") %>% 
                     dplyr::mutate(Readout = "Sequencing")) %>% 
  dplyr::bind_rows(data.table::fread("~/Downloads/INTERNAL_Oncology Reference Screen Info for DMC Partners - AIR.csv") %>% 
                     dplyr::mutate(Readout = "Sequencing-AIR"))




# AIR
df.air <- data.table::fread("~/Downloads/air/filtered_counts_APS005.csv") %>% 
  dplyr::bind_rows(data.table::fread("~/Downloads/air/filtered_counts_APS006.csv")) %>%
  dplyr::bind_rows(data.table::fread("~/Downloads/air/filtered_counts_AIR001.csv")) %>%
  dplyr::rename(CompoundPlate = pert_plate, SampleID = pert_id) 

# APS
df.aps <- data.table::fread("OncRef data/aps/filtered_counts_APS005.csv") %>%
  dplyr::bind_rows(data.table::fread("OncRef data/aps/filtered_counts_APS006.csv")) %>%
  dplyr::rename(CompoundPlate = pert_plate, SampleID = pert_id) 


# MTS
df.mts <- data.table::fread("OncRef data/mts/filtered_counts_MTS28_PMTS086.csv") %>%
  dplyr::bind_rows(data.table::fread("OncRef data/mts/filtered_counts_MTS29_PMTS088.csv")) %>%
  dplyr::bind_rows(data.table::fread("OncRef data/mts/filtered_counts_MTS30_PMTS090.csv")) %>%
  dplyr::bind_rows(data.table::fread("OncRef data/mts/filtered_counts_MTS031.csv")) %>%
  dplyr::rename(CompoundPlate = pert_plate, SampleID = pert_id) 

# CPS
df.cps <- data.table::fread("OncRef data/cps/filtered_counts_cps013.csv") %>%
  dplyr::bind_rows(data.table::fread("OncRef data/cps/filtered_counts_cps014.csv")) %>%
  dplyr::rename(CompoundPlate = pert_plate, SampleID = pert_id) 


# OLD
df.old <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqFilteredCounts.csv") 


df <- df.old %>% 
  dplyr::bind_rows(df.mts) %>% 
  dplyr::bind_rows(df.aps) %>% 
  dplyr::bind_rows(df.cps) %>% 
  dplyr::bind_rows(df.air)



df <- df.old %>% 
  df.mts %>% 
  dplyr::bind_rows(df.aps) %>% 
  dplyr::bind_rows(df.cps) %>% 
  dplyr::bind_rows(df.air)


rm(df.old, df.mts, df.aps, df.air, df.cps)



df %<>%
  dplyr::filter(CompoundPlate %in% dplyr::filter(comprehensive_list, Release != "DoNotRelease")$CompoundPlate) 


df.control <- df %>% 
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon"))


df.trt <- df %>%
  dplyr::filter(!pert_type %in% c("ctl_vehicle", "trt_poscon")) %>% 
  dplyr::mutate(pert_type == "trt_cp") %>%
  dplyr::inner_join(comprehensive_list %>%
                      dplyr::filter(Release != "DoNotRelease") %>% 
                      dplyr::distinct(CompoundPlate, SampleID)) 


df <- df.trt %>%
  dplyr::bind_rows(df.control)


rm(df.trt, df.control)


df %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_type) %>%
  dplyr::mutate(has.data = T) %>% 
  dplyr::full_join(comprehensive_list) %>% 
  dplyr::mutate(has.data = !is.na(has.data)) %>% 
  #dplyr::count(Release, Screen, has.data, Readout) %>% 
  dplyr::filter(has.data, is.na(Release)) %>% 
  View


df.old %>% 
   colnames() %>% 
  paste0(collapse = ", ")
  head


df %<>% 
  dplyr::select(pcr_plate, pcr_well, pert_dose, pert_dose_unit, CompoundPlate, SampleID, pert_type, day, cell_set, bio_rep, replicate_plate, cb_ladder, pert_vehicle, lua, cb_name, cb_log2_dose, n, pool_id, depmap_id, growth_condition) %>% 
  dplyr::distinct()


df %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqFilteredCounts26Q3.csv")




compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundList.csv")


cl <- comprehensive_list %>% 
  dplyr::filter(Readout != "Luminex", Release != "DoNotRelease") %>% 
  dplyr::mutate(Readout == "Sequence") %>% 
  dplyr::distinct(CompoundPlate, SampleID, CompoundName, TargetOrMechanism, ChEMBLID, PubChemCID, GeneSymbolOfTargets, Synonyms,
                  Class, Modality, TargetPathway, HighestClinicalStage, Readout, Release, Prioritized) %>% 
  dplyr::left_join(df %>% 
                     dplyr::distinct(CompoundPlate, SampleID) %>%
                     dplyr::mutate(hasData = TRUE)) %>% 
  dplyr::mutate(hasData = !is.na(hasData),
                Release = word(Release,1,2)) 
  
cl %<>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::mutate(Prioritized = ifelse(n() == 1, TRUE, Prioritized)) %>% 
  dplyr::ungroup() 
                

cl %>% 
  write_csv("OncRef data/PRISMOncologyReferenceSeqCompoundList26Q3.csv")



# Needs to be reviewed! 
cl %>% 
  dplyr::group_by(CompoundName) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::ungroup() %>% 
  View


