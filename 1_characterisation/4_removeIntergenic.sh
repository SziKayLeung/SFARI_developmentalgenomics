LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))

library("data.table")
library("dplyr")

input.file <- "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt"
input.file <- SQANTI_class_preparation(input.file,"nstandard")

input.file.nointergenic <- input.file %>% filter(structural_category != "Intergenic")
write.table(input.file.nointergenic, "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noIntergenic.txt", quote = F, sep = "\t")

unique(input.file.nointergenic$structural_category)
