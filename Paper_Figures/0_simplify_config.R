#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("vroom"))

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"target_gene_annotation/summarise_gene_stats.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)


## ------------ directory names --------------- 

root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
root_rb_dir <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/"
recovered_dir <- "/lustre/recovered/Research_Project-MRC148213/sl693/RBFetal/"
dirnames <- list(
  
  wholetarg_SQ = paste0(root_rb_dir,"6_sqanti3/"),
  
  output = paste0(root_sfari,"/0_output/"),
  utils = paste0(root_sfari,"/0_utils/"),
  protein = paste0(root_sfari, "/8_longReadProteogenomics/longReadProteogenomics"),
  
  DGE = paste0(root_sfari, "10_deseq/1_DGE/"), 
  DTE = paste0(root_sfari, "10_deseq/2_DTE/"),
  DIU = paste0(root_sfari, "10_deseq/3_DIU/")
)


TargetGene = read.table(paste0(root_sfari, "0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]


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
message("Number of transcripts in whole transcriptome dataset after SQANTI filtering, 2 reads 2 samples: ", nrow(class.files$glob_SQ))
message("Number of genes: ", length(unique(class.files$glob_SQ$associated_gene)))

# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons != "Intergenic_1"))
message("Number of transcripts in whole transcriptome dataset after SQANTI filtering, 2 reads 2 samples, -monoexonic intergenic: ", nrow(class.files$glob_SQ))
message("Number of genes: ", length(unique(class.files$glob_SQ$associated_gene)))

# remove mono-exonic transcripts 
refExonNum <- read.csv(paste0(dirnames$utils,"gencode.v40.annotation.numExon.csv"))
exonicClass <- refExonNum %>% 
  group_by(Gene, GeneName, GeneType) %>% 
  dplyr::summarise(maxNumExon = max(NumExon)) %>% 
  mutate(monoExonic = ifelse(maxNumExon == 1, TRUE, FALSE))

refGeneType <- read.table(paste0(dirnames$utils,"gencode.v38.annotation.geneannotation.txt"), header = T)

# identify multi-exonic genes
multiExonicGenes <- unique(exonicClass[exonicClass$monoExonic == FALSE,][["GeneName"]])
# TRUE = monotranscript within multi-exonic gene
class.files <- lapply(class.files, function(x) x %>% mutate(monomulti = ifelse(exons == 1 & associated_gene %in% multiExonicGenes, TRUE,FALSE)))
# filter FALSE to retain multi-transcripts within multi-exonic gene, and mono-transcripts within mono-exonic gene
class.files <- lapply(class.files, function(x) x %>% filter(monomulti == FALSE))
message("Number of transcripts in whole transcriptome dataset after SQANTI filtering, 2 reads 2 samples, -monoexonic intergenic, -monoexonic multigene: ", nrow(class.files$glob_SQ))
message("Number of genes: ", length(unique(class.files$glob_SQ$associated_gene)))

# mmExonicClass <- refExonNum[refExonNum$NumExon == 1 & refExonNum$GeneName %in% exonicClass[exonicClass$monoExonic != TRUE,][["GeneName"]],]
# monomulti.class.files$glob_SQ[monomulti.class.files$glob_SQ$associated_gene %in% mmExonicClass$GeneName,]

# filter targeted dataset to just the target genes
class.files$targ_SQ <- class.files$targ_SQ %>% filter(associated_gene %in% TargetGene)

save(class.files, file = paste0(dirnames$wholetarg_SQ,"all_filtered_classification_2reads2samples_noMonoIntergenicAll.RData"))
write.table(class.files$glob_SQ, paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples_monomultirem.txt"),quote=F, sep = "\t")
write.table(class.files$glob_targ_SQ, paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_2reads2samples_monomultirem.txt"),quote=F, sep = "\t", row.names= F)

## -------------- DESeq2

# Whole DESeq2
WholeDESeqSig <- list(
  sex = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_sex_resSig.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_development_resSig.csv"),delim = ",",show_col_types = FALSE)
)

WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_prenatalsex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_group_removinglateprenatal.csv")))
)

# Expression
Exp <- list(
  whole_group = vroom(paste0(dirnames$DTE,"DESeq2_whole_development_normSig.csv"),delim = ","),
  ns1 = read.csv(paste0(dirnames$DTE,"ONT15_1709_8237_whole_normAll.csv"), header = F)
)
colnames(Exp$ns1) <- colnames(Exp$whole_group)
Exp$whole_group <- rbind(Exp$whole_group,Exp$ns1)
Exp <- lapply(Exp, function(x) merge(x, phenotype, by="sample"))
save(Exp, file = paste0(dirnames$DTE,"DESeq2_whole_normSig.RData"))

ExpGenes <- list(
  whole_sex = fread(paste0(dirnames$DGE,"DESeq2_whole_sex_normAll.csv")),
  whole_group = fread(paste0(dirnames$DGE,"DESeq2_whole_development_normAll.csv"))
)
ExpGenes <- lapply(ExpGenes, function(x) merge(x,phenotype,by="sample",all.x=T))
save(ExpGenes, file = paste0(dirnames$DGE,"DESeq2_whole_normAll.RData"))

ExpGenesSig <- list(
  targeted_sex =  ExpGenes$targeted_sex %>% filter(associated_gene %in% TargetedDESeqGeneSig$sex$associated_gene),
  targeted_group = ExpGenes$targeted_group %>% filter(associated_gene %in% TargetedDESeqGeneSig$age$associated_gene),
  whole_sex = ExpGenes$whole_sex %>% filter(associated_gene %in% WholeDESeqGeneSig$sex$associated_gene),
  whole_group = ExpGenes$whole_group %>% filter(associated_gene %in% WholeDESeqGeneSig$group$associated_gene)
)
save(ExpGenesSig, file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_gene_sig_removinglateprenatal.RData"))