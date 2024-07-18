#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: remove mono-exonic intergenic isoforms
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## --------------------------------


LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))

# packages
suppressMessages({
  library("data.table")
  library("dplyr")
})


# Whole+Targeted merged dataset
sqantiDir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/"
input.file <- paste0(sqantiDir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt")
input.file <- SQANTI_class_preparation(input.file,"nstandard")
input.file.monointergenic <- input.file %>% filter(structural_category == "Intergenic" & exons == 1)
input.file.nomonointergenic <- input.file %>% filter(!isoform %in% input.file.monointergenic$isoform )
write.table(input.file.nomonointergenic, paste0(sqantiDir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt"), quote = F, sep = "\t")

# whole dataset
input.file <- paste0(sqantiDir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples.txt")
input.file <- SQANTI_class_preparation(input.file,"nstandard")
input.file.monointergenic <- input.file %>% filter(structural_category == "Intergenic" & exons == 1)
input.file.nomonointergenic <- input.file %>% filter(!isoform %in% input.file.monointergenic$isoform )
write.table(input.file.nomonointergenic, paste0(sqantiDir,"/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic.txt"), quote = F, sep = "\t")