#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------


LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))


## ------------ directory names --------------- 
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  #glob_root = paste0(root_dir, "rTg4510/A_IsoSeq_Whole"),
  
  # targeted sequencing (Iso-Seq, ONT)
  targ_SQ = paste0(root_dir, "RBFetal/1_SQANTI3Filtered/")

)


## ------------- Phenotype files -------------------



## -------------- Final classification files ------------- 
class.names.files <- list(
  # targeted SQANTI dataset (filtered using relaxed JSON file)
  targ_SQ = paste0(dirnames$targ_SQ, "merged_RulesFilter_result_classification.targetgenes_counts.txt"),
  
  # targeted SQANTI dataset further filtered by minimum 2 reads and 2 samples
  targ_SQ_fil = paste0(dirnames$targ_SQ, "merged_RulesFilter_result_classification.targetgenes_counts_filtered.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))