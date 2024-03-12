#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("vroom"))

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
  
  output = paste0(root_sfari,"0_output/")
)

TargetGene = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/vc362/fetalBrain/genes.txt")[["V1"]]


## ------------- Demultiplex files -------------------

demux.names.files <- list(
  # targeted + whole SQANTI dataset 
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"WholeTargeted_demux_2reads2samples_SQANTIfiltered.csv")
)

demux.files <- lapply(demux.names.files, function(x) fread(x))


## ------------- Phenotype files -------------------
phenotype <- list(
  WholeTargeted = fread(paste0(root_dir, "RBFetal/00_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
)

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
mono.class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons == "Intergenic_1"))
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
TargetedDESeq <- list(
  sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_sex_removinglateprenatal.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_group_removinglateprenatal.csv"),delim = ",",show_col_types = FALSE),
  prenatal_sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_prenatalsex_removinglateprenatal.csv"),delim = ",",show_col_types = FALSE)
)
TargetedDESeq <- lapply(TargetedDESeq, function(x) merge(x, class.files$glob_targ_SQ[,c("isoform","chrom")], by = "isoform"))
save(TargetedDESeq, file = paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_removinglateprenatal.RData"))

TargetedDESeqSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_sex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_group_removinglateprenatal.csv")))
)
TargetedDESeqSig <- lapply(TargetedDESeqSig, function(x) merge(x, class.files$glob_targ_SQ[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory")], by = "isoform", all.x = T))
save(TargetedDESeqSig, file = paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_removinglateprenatal_sig.RData"))

# Whole DESeq2
WholeDESeq <- list(
  sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/anno_whole_sex_removinglateprenatal_TEST.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(root_dir,"RBFetal/4_deseq2/anno_whole_group_removinglateprenatal_TEST.csv"),delim = ",",show_col_types = FALSE)
)
WholeDESeq <- lapply(WholeDESeq, function(x) merge(x, class.files$glob_targ_SQ[,c("isoform","chrom")], by = "isoform"))
save(WholeDESeq, file = paste0(root_dir,"RBFetal/4_deseq2/anno_whole_removinglateprenatal_TEST.RData"))

WholeDESeqSig <- list(
  sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_sex_removinglateprenatal_TEST.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_group_removinglateprenatal_TEST.csv"),delim = ",",show_col_types = FALSE)
)
WholeDESeqSig  <- lapply(WholeDESeqSig , function(x) merge(x, class.files$glob_targ_SQ[,c("isoform","associated_gene","associated_transcript","chrom","structural_category","subcategory")], by = "isoform"))
save(WholeDESeqSig, file = paste0(root_dir,"RBFetal/4_deseq2/anno_whole_removinglateprenatal_TEST_sig.RData"))


WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_prenatalsex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_group_removinglateprenatal.csv")))
)

TargetedDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_gene_prenatalsex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_gene_group_removinglateprenatal.csv")))
)


# Expression 
Exp <- list(
  targeted_sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_targeted_sex_removinglateprenatal.csv"),delim = ","),
  targeted_group = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_targeted_group_removinglateprenatal.csv"),delim = ","),
  whole_sex = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_whole_sex_removinglateprenatal_TEST.csv"),delim = ","),
  whole_group = vroom(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_whole_group_removinglateprenatal_TEST.csv"),delim = ",")
)
Exp <- lapply(Exp, function(x) merge(x, phenotype$WholeTargeted, by="sample"))
Exp <- lapply(Exp, function(x) merge(x, class.files$glob_targ_SQ[,c("isoform","associated_gene","structural_category")], by="isoform", all.x = T))
save(Exp, file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_removinglateprenatal.RData"))

ExpGenes <- list(
  targeted_sex = fread(paste0(root_dir,"RBFetal/4_deseq2/output_norm_targeted_gene_sex_removinglateprenatal.csv")),
  targeted_group = fread(paste0(root_dir,"RBFetal/4_deseq2/output_norm_targeted_gene_group_removinglateprenatal.csv")),
  whole_sex = fread(paste0(root_dir,"RBFetal/4_deseq2/output_norm_whole_gene_sex_removinglateprenatal.csv")),
  whole_group = fread(paste0(root_dir,"RBFetal/4_deseq2/output_norm_whole_gene_group_removinglateprenatal.csv"))
)
ExpGenes <- lapply(ExpGenes, function(x) merge(x,phenotype$WholeTargeted,by="sample",all.x=T))
save(ExpGenes, file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_gene_removinglateprenatal.RData"))

ExpGenesSig <- list(
  targeted_sex =  ExpGenes$targeted_sex %>% filter(associated_gene %in% TargetedDESeqGeneSig$sex$associated_gene),
  targeted_group = ExpGenes$targeted_group %>% filter(associated_gene %in% TargetedDESeqGeneSig$age$associated_gene),
  whole_sex = ExpGenes$whole_sex %>% filter(associated_gene %in% WholeDESeqGeneSig$sex$associated_gene),
  whole_group = ExpGenes$whole_group %>% filter(associated_gene %in% WholeDESeqGeneSig$group$associated_gene)
)
save(ExpGenesSig, file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_gene_sig_removinglateprenatal.RData"))
