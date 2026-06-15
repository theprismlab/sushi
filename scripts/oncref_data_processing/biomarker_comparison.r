library(useful)
library(magrittr)
library(tidyverse)
library(ggrepel)
library(tidytext)
library(scales)
library(ggthemes)


file = "~/code/depmap_data/depmap_24Q4_internal.h5"
source("scripts/biomarker/biomarker_functions.R")


bm <- data.table::fread("OncRef data/bm_lfc_comparison_all.csv")

bm.air <- data.table::fread("OncRef data/bm_lfc_comparison_air_filtered.csv")



bm.air %<>% 
  dplyr::distinct(CompoundName, SampleID) %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::top_n(1, CompoundName) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(bm.air %>% 
                     dplyr::select(-CompoundName, -GeneSymbolOfTargets, -is.target) %>% 
                     dplyr::distinct())
  

bm %<>% 
  dplyr::distinct(CompoundName, SampleID) %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::top_n(1, CompoundName) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(bm %>% 
                     dplyr::select(-CompoundName, -GeneSymbolOfTargets, -is.target) %>% 
                     dplyr::distinct())


compound_annotations <- data.table::fread("OncRef data/PRISMOncologyReferenceSeqCompoundListExtended.csv")

targets <- compound_annotations %>%
  dplyr::distinct(SampleID, GeneSymbolOfTargets) %>% 
  separate_rows(GeneSymbolOfTargets, sep = fixed(";")) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(GeneSymbolOfTargets != "N/A") %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::summarise(GeneSymbolOfTargets = paste0(sort(unique(GeneSymbolOfTargets)), collapse = ";")) %>% 
  dplyr::ungroup()
  

tr <- bm %>% 
  dplyr::distinct(SampleID, feature) %>%
  dplyr::left_join(targets) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(is.target = (word(feature, sep = fixed("--")) == GeneSymbolOfTargets) | 
                  (word(feature, sep = fixed(".")) == GeneSymbolOfTargets)) %>% 
  dplyr::group_by(SampleID, feature) %>% 
  dplyr::summarise(is.target = any(is.target, na.rm = T)) %>% 
  dplyr::ungroup()

bm %<>% 
  dplyr::left_join(tr)



tr <- bm.air %>% 
  dplyr::distinct(SampleID, feature) %>%
  dplyr::left_join(targets) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(is.target = (word(feature, sep = fixed("--")) == GeneSymbolOfTargets) | 
                  (word(feature, sep = fixed(".")) == GeneSymbolOfTargets)) %>% 
  dplyr::group_by(SampleID, feature) %>% 
  dplyr::summarise(is.target = any(is.target, na.rm = T)) %>% 
  dplyr::ungroup()

bm.air %<>% 
  dplyr::left_join(tr)



bm.air %>%
  dplyr::filter(feature_set == "Expression", 
                type == "no.SR", correlation_coef < 0,
                response %in% c("corrected", "uncorrected")) %>%
  dplyr::group_by(cn, response) %>% 
  dplyr::arrange(desc(rank)) %>% dplyr::mutate(dc = -c(0, diff(correlation_coef))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.target, rank == 1) %>%
  dplyr::rename(cor = correlation_coef) %>% 
  dplyr::select(-p_val, -q_val, -n, -type, -is.target) %>% 
  tidyr::pivot_wider(names_from = response,
                     values_from = c("rank", "cor", "dc")) %>% 
  dplyr::filter(is.finite(rank_corrected), is.finite(rank_uncorrected)) %>% 
  ggplot() +
  geom_point(aes(x = abs(cor_corrected), y = abs(cor_uncorrected))) +
  geom_abline() +
  facet_wrap(feature + CompoundName ~ .)


bm.air %>%
  dplyr::filter(feature_set == "Expression", 
                type == "no.SR", correlation_coef < 0,
                response %in% c("corrected", "uncorrected")) %>%
  dplyr::group_by(cn, response) %>% 
  dplyr::arrange(desc(rank)) %>% dplyr::mutate(dc = -c(0, diff(correlation_coef))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.target, rank == 1) %>%
  dplyr::rename(cor = correlation_coef) %>% 
  dplyr::select(-p_val, -q_val, -n, -type, -is.target) %>% 
  tidyr::pivot_wider(names_from = response,
                     values_from = c("rank", "cor", "dc")) %>% 
  dplyr::filter(is.finite(rank_corrected), is.finite(rank_uncorrected)) %>% 
  ggplot() +
  geom_point(aes(x = dc_corrected, y = dc_uncorrected)) +
  geom_abline() +
  facet_wrap(feature + CompoundName ~ .)
  


# AIR FILTERING COMPARISON ----


bm.air %>% 
  dplyr::filter(!is.target, rank <= 10) %>% 
  View

bm.counts <- bm.air %>% 
  dplyr::filter(!is.target, rank <= 10) %>% 
  dplyr::group_by(type, response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(type, response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup()

  
bm.counts %>% 
  dplyr::distinct(type, response, feature_set, feature, n.c) %>%
  dplyr::group_by(type, response) %>% 
  dplyr::top_n(10, n.c) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(bm.counts) %>%
  dplyr::filter(response %in% c("corrected", "uncorrected")) %>% 
  dplyr::mutate(response = ifelse(response == "corrected", "Corrected", " Uncorrected"), 
                type = ifelse(type == "full", "All",
                              ifelse(type == "no.S", "AS removed", "AS/R removed"))) %>% 
  ggplot() +
  geom_bar(aes(x = reorder_within(paste0(substr(feature_set,1,3), "_", feature), -n.c, paste0(type, response)), 
               fill = as.factor(n))) +
  facet_wrap(response ~ type, scales = "free_x", ncol = 3) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=70, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\nDoses")
  


bm.counts %>% 
  dplyr::distinct(type, response, feature_set, feature, n.c) %>%
  dplyr::group_by(type, response) %>% 
  dplyr::top_n(15, n.c) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(bm.counts) %>%
  dplyr::filter(type %in% c("no.S", "no.SR")) %>% 
  dplyr::mutate(type = ifelse(type == "full", "All",
                              ifelse(type == "no.S", "AS removed", "AS/R removed"))) %>% 
  ggplot() +
  geom_bar(aes(x = reorder_within(paste0(substr(feature_set,1,2), "_", feature), -n.c, paste0(type, response)), 
               fill = as.factor(n))) +
  facet_wrap(type ~ response, scales = "free_x", nrow = 2) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=70, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\nDoses")



bm.air %>%
  dplyr::filter(is.target,
                feature_set == "Expression",
                rank <= 1,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(type, response, CompoundPlate, CompoundName, feature_set, feature) %>%
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(type) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(type, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps.long = paste(sort(unique(CompoundName)), collapse = "\n"),
                   cmps = paste(sort(word(unique(CompoundName), sep = fixed("-"))), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = cmps), # ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(paste0(type, " (", n.c,  ")") ~ ., ncol = 3) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)")

bm.air %>%
  dplyr::filter(is.target,
                feature_set == "Expression",
                rank <= 10,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(type, response, CompoundPlate, CompoundName, feature_set, feature) %>%
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(type) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(type, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps.long = paste(sort(unique(CompoundName)), collapse = "\n"),
                   cmps = paste(sort(word(unique(CompoundName), sep = fixed("-"))), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = cmps), # ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(paste0(type, " (", n.c,  ")") ~ ., ncol = 3) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)")

bm.air %>%
  dplyr::filter(is.target,
                feature_set == "Expression",
                rank <= 100,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(type, response, CompoundPlate, CompoundName, feature_set, feature) %>%
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(type) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(type, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps.long = paste(sort(unique(CompoundName)), collapse = "\n"),
                   cmps = paste(sort(word(unique(CompoundName), sep = fixed("-"))), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = cmps), # ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(paste0(type, " (", n.c,  ")") ~ ., ncol = 3) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)")
  

# Let's extend the set to larger set ----

bm %<>% 
  dplyr::filter(!SampleID %in% bm.air$SampleID)
  

bm.counts <- bm %>% 
  dplyr::filter(response != "adaptive_mild") %>% 
  dplyr::filter(!is.target, rank <= 1) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup()


p1 = bm.counts %>% 
  dplyr::distinct(response, feature_set, feature, n.c) %>%
  dplyr::group_by(response) %>% 
  dplyr::top_n(10, n.c) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(bm.counts) %>%
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_mild", "adaptive_heavy")) %>% 
  ggplot() +
  geom_bar(aes(x = reorder_within(paste0(substr(feature_set,1,3), "_", feature), -n.c, response), 
               fill = as.factor(n))) +
  facet_wrap(word(response, sep = fixed("_")) ~ ., scales = "free_x", nrow = 1) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=70, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\nDoses", title = "Common Top 1 Non-target Correlations")

bm.counts <- bm %>% 
  dplyr::filter(!is.target, rank <= 10) %>% 
  dplyr::filter(response != "adaptive_mild") %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup()


p2 = bm.counts %>% 
  dplyr::distinct(response, feature_set, feature, n.c) %>%
  dplyr::group_by(response) %>% 
  dplyr::top_n(10, n.c) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(bm.counts) %>%
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_mild", "adaptive_heavy")) %>% 
  ggplot() +
  geom_bar(aes(x = reorder_within(paste0(substr(feature_set,1,3), "_", feature), -n.c, response), 
               fill = as.factor(n))) +
  facet_wrap(word(response, sep = fixed("_")) ~ ., scales = "free_x", nrow = 1) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=70, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\nDoses", title = "Common Top 10 Non-target Correlations")



cowplot::plot_grid(p1,p2, ncol = 1)

bm %>%
  dplyr::filter(is.target,
                #feature_set %in% c("Expression"), # "CRISPR", "RNAi"),
                rank <= 1,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(response, CompoundPlate, CompoundName, feature_set, feature) %>% 
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(feature_set) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(feature_set, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps = paste(sort(unique(CompoundName)), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(  . ~ paste0(feature_set, " (", n.c,  ")")) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)",
       title = "Top 1 Target Recovery")



bm %>%
  dplyr::filter(is.target,
                #feature_set %in% c("Expression"), # "CRISPR", "RNAi"),
                rank <= 10,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(response, CompoundPlate, CompoundName, feature_set, feature) %>% 
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(feature_set) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(feature_set, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps = paste(sort(unique(CompoundName)), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(  . ~ paste0(feature_set, " (", n.c,  ")")) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)",
       title = "Top 10 Target Recovery")



bm %>%
  dplyr::filter(is.target,
                #feature_set %in% c("Expression"), # "CRISPR", "RNAi"),
                rank <= 100,
                response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(response, CompoundPlate, CompoundName, feature_set, feature) %>% 
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>% 
  dplyr::group_by(feature_set) %>% dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::group_by(feature_set, corrected, uncorrected, n.c) %>% 
  dplyr::summarise(n = length(unique(CompoundName)), 
                   cmps = paste(sort(unique(CompoundName)), collapse = "\n")) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             color = as.factor(sign(corrected - uncorrected)))) +
  geom_abline(color = "grey") + 
  geom_text(aes( label = n), show.legend = F) +
  geom_text_repel(aes(label = ifelse(abs(corrected - uncorrected) > 0, cmps, NA)),
                  size = 2.5, show.legend = F) + 
  scale_color_wsj() + 
  facet_wrap(  . ~ paste0(feature_set, " (", n.c,  ")")) +
  theme_bw() +
  labs(x = "Number of doses (Corrected)", y = "Number of doses (Uncorrected)",
       title = "Top 100 Target Recovery")




bm %>% # bm.air %>% dplyr::filter(type == "no.S") %>% 
  dplyr::filter(is.target, rank <= 1) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) > 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 1 Multiple Compounds")


bm %>% # bm.air %>% dplyr::filter(type == "no.S") %>% 
  dplyr::filter(is.target, rank <= 10) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) > 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 10 Multiple Compounds")

bm %>% # bm.air %>% dplyr::filter(type == "no.S") %>% 
  dplyr::filter(is.target, rank <= 100) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) > 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 100 Multiple Compounds")


bm %>% 
  dplyr::filter(is.target, rank <= 1) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) == 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(reorder(paste0(feature_set, "_", feature), as.numeric(as.factor(response)), min) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", , fill = "# of\n Doses", title = "Top 1 Single Compounds")



bm %>% 
  dplyr::filter(is.target, rank <= 10) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) == 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(reorder(paste0(feature_set, "_", feature), as.numeric(as.factor(response)), min) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", , fill = "# of\n Doses", title = "Top 10 Single Compounds")



bm %>% 
  dplyr::filter(is.target, rank <= 100) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) == 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(reorder(paste0(feature_set, "_", feature), as.numeric(as.factor(response)), min) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", , fill = "# of\n Doses", title = "Top 100 Single Compounds")




bm.air %>% dplyr::filter(#type == "no.SR", 
                              feature_set == "Expression") %>% 
  dplyr::filter(is.target, rank <= 1) %>% 
  dplyr::group_by(type, response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(type, response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(type, feature_set, feature) %>% 
  dplyr::mutate( flag = max(n.c) - min(n.c) > 0 | length(unique(response)) < 2) %>% 
  dplyr::group_by(type, response, feature_set, feature) %>% 
  dplyr::mutate( cmp =  paste0(unique(word(substr(CompoundName, 1,6), sep = fixed("-"))), collapse = "\n")) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(any(flag)) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  geom_label(aes(x = response, y = 1, label = cmp), alpha = 0.5, size = 3) + 
  facet_grid(type ~ feature) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 1")


p1

p2 = bm.air %>% dplyr::filter(type == "no.SR", 
                              feature_set == "Expression") %>% 
  dplyr::filter(is.target, rank <= 10) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 2) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 10")



p3 = bm.air %>% dplyr::filter(type == "no.SR",
                              feature_set == "Expression") %>% 
  dplyr::filter(is.target, rank <= 100) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 2) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 100")



p1 

cowplot::plot_grid(p1,p2,p3, nrow = 1)


bm.air %>% dplyr::filter(type == "no.SR") %>% 
  dplyr::filter(is.target, rank <= 10) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(paste0(feature_set, "_", feature) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", fill = "# of\n Doses", title = "Top 10")


bm %>% 
  dplyr::filter(is.target, rank <= 1) %>% 
  dplyr::group_by(response, feature_set, feature, CompoundPlate, CompoundName, SampleID) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(response, feature_set, feature) %>%
  dplyr::mutate(n.c = length(unique(CompoundName))) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy")) %>% 
  dplyr::mutate(response = word(response, sep = fixed("_"))) %>% 
  dplyr::group_by(feature_set, feature) %>% 
  dplyr::filter(max(n.c) - min(n.c) > 0 | length(unique(response)) < 3,
                max(n.c) == 1) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = response, 
               fill = as.factor(n))) +
  facet_wrap(reorder(paste0(feature_set, "_", feature), as.numeric(as.factor(response)), min) ~ .) +
  scale_x_reordered() +
  theme_bw() + scale_fill_ptol() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        legend.position = "right") +
  labs(x = NULL, y = "Number of agents", , fill = "# of\n Doses", title = "Top 1 Single Compounds")



# -----

bm = bm0






bm.counts <- bm %>%
#   dplyr::filter(substr(CompoundPlate, 1,4) == "PAPS") %>% 
  dplyr::filter(!is.target, rank <= 1) %>% 
  dplyr::group_by(response, feature_set, feature) %>% 
  dplyr::summarise(n.c = length(unique(CompoundName)),
                   n.t = length(unique(GeneSymbolOfTargets)),
                   n.s = length(unique(cn))) 


bm.counts %>% 
  dplyr::group_by(response, feature) %>% dplyr::mutate(n.cc = sum(n.c)) %>% dplyr::ungroup() %>% 
  dplyr::group_by(response) %>%
  dplyr::top_n(20, n.cc) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(response %in% c("corrected", "uncorrected", "adaptive_heavy", "adaptive_mild")) %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(feature, -n.cc, response), y = n.c, fill = feature_set)) +
  facet_wrap(response ~ ., scales = "free_x",nrow = 1) +
  scale_x_reordered() +
  theme_bw() + scale_fill_wsj() +
  theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1),
        legend.position = "bottom") +
  labs(x = NULL, y = "Number of compounds", fill = NULL)


temp <- bm %>% 
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, CompoundPlate, feature_set,feature, response, rank) %>%
  dplyr::filter(rank <= 10, response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(CompoundPlate, CompoundName, response, feature_set, feature) %>% 
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) 




temp2 <- bm %>% 
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, CompoundPlate, feature_set,feature, response, rank) %>%
  dplyr::filter(rank <= 10, response %in% c("corrected", "uncorrected")) %>% 
  dplyr::count(CompoundPlate, CompoundName, response, feature_set, feature) %>%
  dplyr::group_by(CompoundPlate, CompoundName, response) %>% 
  dplyr::top_n(1, n) %>% 
  dplyr::distinct(CompoundPlate, CompoundName, response, n) %>%
  tidyr::pivot_wider(names_from = response, values_from = n, values_fill = 0) %>%
  dplyr::ungroup()





temp %>%
  # dplyr::filter(grepl("AIR", CompoundName)) %>% View
  dplyr::count(feature_set, corrected, uncorrected) %>% 
  ggplot(aes(x = as.factor(corrected), y = as.factor(uncorrected), label = n, fill = n)) +
  geom_tile(show.legend = F) + 
  geom_text(show.legend = F) +
  facet_wrap(feature_set ~ .) +
  theme_few() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Number of doses (corrected)", y = "Number of doses (uncorrected)",
       title = "Compound x Target pairs rank in top 10")


temp %>% 
  dplyr::filter(abs(uncorrected - corrected) > 1) %>% 
  dplyr::distinct(CompoundName) %>% 
  dplyr::left_join(temp) %>% 
  ggplot() +
  geom_point(aes(y = uncorrected, x = corrected, color = feature_set)) +
  geom_abline() + 
  facet_wrap(paste0(CompoundName) ~ .)
   head





temp2 %>% 
  dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::count(corrected, uncorrected) %>%
  ggplot(aes(x = as.factor(corrected), y = as.factor(uncorrected), label = n, fill = n)) +
  geom_tile(show.legend = F) + 
  geom_text(show.legend = F) +
  theme_few() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Number of doses (corrected)", y = "Number of doses (uncorrected)",
       title = "Compound x Target Pairs")









temp %>%
  dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::count(feature_set, corrected, uncorrected) %>% 
  ggplot(aes(x = as.factor(corrected), y = as.factor(uncorrected), label = n, fill = n)) +
  geom_tile(show.legend = F) + 
  geom_text(show.legend = F) +
  facet_wrap(feature_set ~ .) +
  theme_few() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Number of doses (corrected)", y = "Number of doses (uncorrected)",
       title = "Compound x Target pairs rank in top 10 (AiR only)")

# We have 65 air agents, 49 ranks at least one target in top 100 ranks, 46 in top 10, 36 in top 1.

# Corrected:     We have 65 air agents, 47 ranks at least one target in top 100 ranks, 45 in top 10, 34 in top 1.
# Uncorrected:   We have 65 air agents, 47 ranks at least one target in top 100 ranks, 42 in top 10, 33 in top 1.
# Adaptive:      We have 65 air agents, 48 ranks at least one target in top 100 ranks, 45 in top 10, 34 in top 1.


# Corrected:     We have 237 air agents, 172 ranks at least one target in top 100 ranks, 159 in top 10, 129 in top 1.
# Uncorrected:   We have 237 air agents, 177 ranks at least one target in top 100 ranks, 155 in top 10, 127 in top 1.
# Adaptive:      We have 237 air agents, 173 ranks at least one target in top 100 ranks, 160 in top 10, 130 in top 1.


bm %>% 
  # dplyr::distinct(CompoundName) %>% dplyr::mutate(is.air = grepl("AIR", CompoundName)) %>%
  # dplyr::filter(is.air) %>%
  # dplyr::left_join(bm) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.target, rank <= 1, response == "corrected") %>% 
  dplyr::count(CompoundName) 
  


temp2 %>% 
  ggplot(aes(x = corrected, y = uncorrected,
             label = ifelse(corrected != uncorrected, CompoundName, NA))) +
  geom_abline() + 
  geom_point() +
  geom_text_repel(max.overlaps = 100, size = 3) +
  theme_bw()
  

dplyr::count(feature_set, corrected, uncorrected) 




bm %>% 
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")),
                bin2 = ifelse( adaptive_mild > 10, 
                               ifelse(adaptive_heavy > 10, "none", "heavy"),
                               ifelse(adaptive_heavy > 10, "mild", "both")),
                bin3 = ifelse( adaptive_mild_complete > 10, 
                               ifelse(adaptive_heavy_complete > 10, "none", "heavy"),
                               ifelse(adaptive_heavy_complete > 10, "mild", "both")))%>% 
  dplyr::group_by(bin, bin2) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>%
  ggplot() +
  geom_abline(lty = 2) + 
  geom_jitter(aes(x = adaptive_mild, 
                  y = adaptive_heavy), shape = 1) +
  geom_text(aes(x = ifelse(bin2 %in% c("both", "mild"), 3, 30),
                y =  ifelse(bin2 %in% c("both", "heavy"), 3, 30), label = nn),
            color = "red") + 
  scale_y_log10(limits = c(1,100)) + scale_x_log10(limits = c(1,100)) +
  facet_wrap(bin ~ ., scales = "free") 


bm %>%
  dplyr::filter(is.target) %>%
  dplyr::rowwise() %>% dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::distinct(cn, CompoundName, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")),
                bin2 = ifelse( adaptive_mild > 10, 
                               ifelse(adaptive_heavy > 10, "none", "heavy"),
                               ifelse(adaptive_heavy > 10, "mild", "both")),
                bin3 = ifelse( adaptive_mild_complete > 10, 
                               ifelse(adaptive_heavy_complete > 10, "none", "heavy"),
                               ifelse(adaptive_heavy_complete > 10, "mild", "both")))%>% 
  dplyr::group_by(bin, bin2) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>%
  ggplot() +
  geom_abline(lty = 2) + 
  geom_jitter(aes(x = adaptive_mild, 
                  y = adaptive_heavy), shape = 1) +
  geom_text(aes(x = ifelse(bin2 %in% c("both", "mild"), 3, 30),
                y =  ifelse(bin2 %in% c("both", "heavy"), 3, 30), label = nn),
            color = "red") + 
  scale_y_log10(limits = c(1,100)) + scale_x_log10(limits = c(1,100)) +
  facet_wrap(bin ~ ., scales = "free") 


bm %>%
  dplyr::filter(is.target) %>%
  dplyr::rowwise() %>% dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::distinct(cn, CompoundName, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")),
                bin2 = ifelse( adaptive_heavy_complete > 10, 
                               ifelse(adaptive_heavy > 10, "none", "heavy"),
                               ifelse(adaptive_heavy > 10, "heavy_complete", "both")))%>% 
  dplyr::group_by(bin, bin2) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>%
  ggplot() +
  geom_abline(lty = 2) + 
  geom_jitter(aes(x = adaptive_heavy_complete, 
                  y = adaptive_heavy), shape = 1) +
  geom_text(aes(x = ifelse(bin2 %in% c("both", "heavy_complete"), 3, 30),
                y =  ifelse(bin2 %in% c("both", "heavy"), 3, 30), label = nn),
            color = "red") + 
  scale_y_log10(limits = c(1,100)) + scale_x_log10(limits = c(1,100)) +
  facet_wrap(bin ~ ., scales = "free") 





bm %>%
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, SampleID, CompoundPlate, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")),
                bin2 = ifelse( adaptive_heavy > 10, "missed", "adaptive"))%>% 
  dplyr::group_by(bin, bin2) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>%
  dplyr::count(feature_set, feature, CompoundName, CompoundPlate, SampleID,  bin, bin2) %>% 
  tidyr::pivot_wider(names_from = bin2, values_from = n, values_fill = 0) %>% 
  dplyr::filter(bin != "none") %>% # , bin != "both") %>% 
  # dplyr::rowwise() %>% dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::group_by(feature_set, feature, CompoundName, CompoundPlate, SampleID) %>% 
  dplyr::summarise(bins = paste0(sort(unique(bin[missed + adaptive == max(missed + adaptive)])), collapse = ";"),  #paste0(sort(unique(bin)), collapse = ";"), #ifelse(all(bin == "corrected"), "corrected", ifelse(all(bin == "uncorrected"), "uncorrected", "corrected\nuncorrected")), 
                   missed = sum(missed), adaptive = sum(adaptive),
                   ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rate = adaptive/(adaptive + missed),
                sum = (adaptive + missed),
                flag = rate < 1 & sum > 1 & missed > 1) %>% 
  ggplot(aes(y =  rate,
             x = sum,
             color = bins)) +
  geom_jitter(aes(x = ifelse(flag, NA, sum)), height = 0.02, width = 0.1) + 
  geom_point(aes(x = ifelse(flag, sum , NA))) +
  geom_text_repel(aes(label = ifelse(flag,
                                     paste0(CompoundName, "\n", feature), NA)), show.legend = F,
                  size = 3, max.overlaps = 100) + 
  facet_wrap(feature_set ~ .) +
  theme_bw() + scale_color_wsj() + 
  labs(x = "Number of recovered doses in corrected or uncorrected data", 
       y = "Fraction of recovered doses in adaptive method")
  



aug <- bm %>%
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, SampleID, CompoundPlate, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")),
                bin2 = ifelse( adaptive_heavy > 10, "missed", "adaptive"))%>% 
  dplyr::group_by(bin, bin2) %>% dplyr::mutate(nn = n()) %>% 
  dplyr::ungroup() %>%
  dplyr::count(feature_set, feature, CompoundName, CompoundPlate, SampleID,  bin, bin2) %>% 
  tidyr::pivot_wider(names_from = bin2, values_from = n, values_fill = 0) %>% 
  dplyr::filter(bin != "none") %>% # , bin != "both") %>% 
  # dplyr::rowwise() %>% dplyr::filter(grepl("AIR", CompoundName)) %>% 
  dplyr::group_by(feature_set, feature, CompoundName, CompoundPlate, SampleID) %>% 
  dplyr::summarise(bins = paste0(sort(unique(bin[missed + adaptive == max(missed + adaptive)])), collapse = ";"),  #paste0(sort(unique(bin)), collapse = ";"), #ifelse(all(bin == "corrected"), "corrected", ifelse(all(bin == "uncorrected"), "uncorrected", "corrected\nuncorrected")), 
                   missed = sum(missed), adaptive = sum(adaptive),
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rate = adaptive/(adaptive + missed),
                sum = (adaptive + missed),
                flag = rate < 1 & sum > 1 & missed > 1) %>% 
  dplyr::distinct(CompoundPlate, SampleID, CompoundName, feature_set, feature, rate, sum, flag)
  
aug %>% 
  dplyr::filter(rate == 0) %>%
  dplyr::distinct(CompoundName, CompoundPlate) %>%
  dplyr::left_join(aug) %>% 
  View
  head

bm %>%
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, CompoundPlate, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::filter(corrected <= 10 | uncorrected <= 10) %>%  
  dplyr::left_join(aug) %>% 
  ggplot() +
  geom_jitter(aes(x = log10(corrected / uncorrected), 
                 y =  log10(adaptive_heavy), #  / pmin(corrected, uncorrected)),
                 color = ifelse(flag, CompoundName, NA))) + 
  theme_bw() + # scale_color_wsj() + 
  facet_wrap(feature_set ~ .) +
  labs(color = NULL) 


bm %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, SampleID, CompoundPlate, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::mutate(bin = ifelse( corrected > 10, 
                              ifelse(uncorrected > 10, "none", "uncorrected"),
                              ifelse(uncorrected > 10, "corrected", "both")))%>% 
  dplyr::count(bin, CompoundName, CompoundPlate, feature_set, feature) %>%
  tidyr::pivot_wider(names_from = bin, values_from = n, values_fill = 0) %>%
  dplyr::arrange(corrected - uncorrected) %>% 
  View
  

bm %>%
  dplyr::filter(is.target) %>% 
  dplyr::distinct(cn, CompoundName, feature_set,feature, response, rank) %>% 
  tidyr::pivot_wider(names_from = response, values_from = rank, values_fill = 100) %>% 
  dplyr::filter(corrected <= 10 | uncorrected <= 10) %>% 
  ggplot() +
  geom_point(aes(x = (corrected / uncorrected), 
                 y = adaptive_heavy)) +
  scale_y_log10() + scale_x_log10()



  
  