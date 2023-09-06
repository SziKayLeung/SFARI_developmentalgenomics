#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------


LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"target_gene_annotation/summarise_gene_stats.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)



## ------------ directory names --------------- 
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  glob_SQ = paste0(root_dir, "/RBFetal/WholeTranscriptome/"),
  
  # targeted sequencing (Iso-Seq, ONT)
  targ_SQ = paste0(root_dir, "RBFetal/1_SQANTI3Filtered/test/")

)


## ------------- Phenotype files -------------------

phenotype <- fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/00_metadata/WholeTargetedphenotype.csv', data.table=F, stringsAsFactors=F)

## -------------- Final classification files ------------- 
class.names.files <- list(
  # whole transcriptome SQANTI dataset
  glob_SQ = paste0(dirnames$glob_SQ, "SQANTI3_whole_RulesFilter_result_classification.counts_filtered.txt"),
  
  # targeted SQANTI dataset (filtered using relaxed JSON file)
  targ_SQ = paste0(dirnames$targ_SQ, "SQANTI3_collapse_options_RulesFilter_result_classification.targetgenes_counts.txt"),
  
  # targeted SQANTI dataset further filtered by minimum 2 reads and 2 samples
  targ_SQ_fil = paste0(dirnames$targ_SQ, "SQANTI3_collapse_options_RulesFilter_result_classification.targetgenes_counts_filtered.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))


## ----------- read lengths ----------------
tpfpostnatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/0_lengths/targeted_postfilter_postnatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
tpfprenatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/targeted_postfilter_prenatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
tpostnatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/targeted_postnatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
tprenatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/targeted_prenatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
wpfpostnatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/whole_postfilter_postnatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
wpfprenatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/whole_postfilter_prenatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
wpostnatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/whole_postnatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)
wprenatal<-fread('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/whole_prenatal_lengths.txt.gz', stringsAsFactors = F, data.table = F)

