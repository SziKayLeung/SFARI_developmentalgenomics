#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))

LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"target_gene_annotation/summarise_gene_stats.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)


## ------------ directory names --------------- 
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/"
root_rb_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  glob_SQ = paste0(root_dir, "/RBFetal/WholeTranscriptome/"),
  
  # targeted sequencing (Iso-Seq, ONT)
  targ_SQ = paste0(root_dir,"RBFetal/1_SQANTI3/"),
  
  wholetarg_SQ = paste0(root_rb_dir,"SQANTI/"),
  
  output = paste0(root_sfari,"/0_output")
)


TargetGene = read.table(paste0(root_sfari, "/0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]
TargetGene20 = read.csv(paste0(root_sfari, "/0_metadata/20SexDifferenceTargetGenes.csv"), header = F)[["V1"]]


## ------------- Phenotype files -------------------
phenotype <- list(
  WholeTargeted = read.csv(paste0(root_sfari, "/12_deseq2/WholeTargetedphenotype.csv")) %>% mutate(time = age)
)

## ------------- Demultiplex files -------------------

demux.names.files <- list(
  # targeted + whole SQANTI dataset 
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"WholeTargeted_demux_2reads2samples_SQANTIfiltered.csv")
)

demux.files <- lapply(demux.names.files, function(x) fread(x))


## -------------- Final classification files ------------- 
class.names.files <- list(
  
  # targeted + whole SQANTI dataset futher filtered by minimum 2 reads and 2 samples
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt"),
  
  # whole transcriptome SQANTI dataset further filtered by minimum 2 reads and 2 samples
  glob_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples.txt"),
  
  # targeted SQANTI dataset further filtered by minimum 2 reads and 2 samples
  targ_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples.txt")
  
) 

class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons != "Intergenic_1"))

# merge class.files with demux for
class.files$glob_targ_SQ <- merge(class.files$glob_targ_SQ, demux.files$glob_targ_SQ, by = "isoform", all.x=T)

# filter targeted dataset to just the target genes
class.files$targ_SQ <- class.files$targ_SQ %>% filter(associated_gene %in% TargetGene)

wholesamples <- colnames(class.files$glob_targ_SQ)[grepl("Whole", colnames(class.files$glob_targ_SQ))]
targetedsamples <- colnames(class.files$glob_targ_SQ)[grepl("Targeted", colnames(class.files$glob_targ_SQ))]
matchedsamples <- intersect(gsub("^.*?Whole","",wholesamples),gsub("^.*?Targeted","",targetedsamples))
targetedmatchedsamples <- paste0("Targeted",matchedsamples)
wholematchedsamples <- paste0("Whole",matchedsamples)

## -------------- DESeq2 
TargetedDESeq2 <- list(
  sex = fread(paste0(root_dir,"/RBFetal/4_deseq2/anno_targeted_sex.csv")),
  age = fread(paste0(root_dir,"/RBFetal/4_deseq2/anno_targeted_group.csv")),
  sex_age = fread(paste0(root_dir,"/RBFetal/4b_deseq2_covariate/anno_targeted_group.csv")),
  sex_age_off = fread(paste0(root_dir,"/RBFetal/4b_deseq2_covariate_off/anno_targeted_group.csv"))
)

WholeDESeq2 <- list(
  sex = fread(paste0(root_sfari,"/12_deseq2/output_whole_sex.csv")),
  age = fread(paste0(root_sfari,"/12_deseq2/output_whole_group.csv")),
)
WholeDESeq2 <- lapply(WholeDESeq2, function(x) merge(x, class.files$glob_targ_SQ, by = "isoform"))


Exp <- list(
  targeted = fread(paste0(root_dir,"RBFetal/4_deseq2/output_norm_targeted_sex.csv")),
  whole = fread(paste0(root_sfari,"/12_deseq2/output_norm_whole_sex.csv"))
)
Exp$targeted <- merge(Exp$targeted, phenotype$WholeTargeted, by="sample")
Exp$targeted <- merge(Exp$targeted, class.files$glob_targ_SQ, by = "isoform", all.x = T)
Exp$wholeanno <- Exp$whole %>% filter(isoform %in% class.files$glob_targ_SQ$isoform)
Exp$wholeanno <- merge(Exp$wholeanno, phenotype$WholeTargeted, by="sample")
Exp$wholeanno <- merge(Exp$wholeanno, class.files$glob_targ_SQ, by = "isoform", all.x = T)

gtf <- list(
  glob_targ = rtracklayer::import(paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))

