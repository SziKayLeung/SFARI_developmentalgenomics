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
LOGEN_ROOT <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"compare_datasets/dataset_identifer.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)


## ------------ directory names --------------- 

root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
dirnames <- list(
  
  # general
  output = paste0(root_sfari,"/0_output/"),
  utils = paste0(root_sfari,"/0_utils/"),

  whole = paste0(root_sfari,"A_Whole/"),
  wholetarg_SQ = paste0(root_sfari,"C_Whole_Targeted/9_sqanti_final/"),
  protein = paste0(root_sfari, "C_Whole_Targeted/10_longReadProteogenomics/"),
  
  # whole dataset differential expression analysis
  DGE = paste0(root_sfari, "A_Whole/10_deseq/1_DGE/"), 
  DTE = paste0(root_sfari, "A_Whole/10_deseq/2_DTE/"),
  DIU = paste0(root_sfari, "A_Whole/10_deseq/3_DIU/")
)


TargetGene = read.table(paste0(root_sfari, "0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]

## ------------- Phenotype files -------------------

phenotype <- fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
phenotype <- phenotype %>% mutate(type = ifelse(grepl("Targeted",sample),"Targeted","Whole"),
                     sampleID = gsub("^Targeted", "", sample)) %>% mutate(sampleID = gsub("^Whole","", sampleID))
phenotype <- phenotype %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), 
                                  col = paste0(sample,"_",group),
                                  sex = factor(sex, levels = c("M","F"))) 
femaleWhole <- phenotype[phenotype$sex == "F" & grepl("Whole",phenotype$sample),][["sample"]]
maleWhole <- phenotype[phenotype$sex == "M" & grepl("Whole",phenotype$sample),][["sample"]]
postWhole <- phenotype[phenotype$group == "Postnatal" & grepl("Whole",phenotype$sample),][["sample"]]
preWhole <- phenotype[phenotype$group == "Prenatal" & grepl("Whole",phenotype$sample),][["sample"]]
matchedsamples <- intersect(phenotype[phenotype$type == "Whole","sampleID"],phenotype[phenotype$type == "Targeted","sampleID"])
wholematchedsamples <- phenotype[phenotype$sampleID %in% matchedsamples & phenotype$type == "Whole","sample"]
targetedmatchedsamples <- phenotype[phenotype$sampleID %in% matchedsamples & phenotype$type == "Targeted","sample"]

# manifest 
manifest <- fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_manifest.csv"),data.table=F, stringsAsFactors=F)

## -------------- Final classification files -------------

## Transcript level: targeted + whole SQANTI dataset file generated from filter 
# keeping only the transcripts kept from sqanti filtering (relaxed json file)
# removed mono-exonic intergenic transcripts 
# removed mono-exonic transcripts within multi-exonic genes
# futher filtered by minimum 2 reads and 2 samples

class.names.files <- list(
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"sqantifiltered_monoexonicfiltered_2reads2samples_classification_finalversion.txt"),
  glob_collapsed = paste0(dirnames$whole,"6_sqanti/Whole_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt")
)
class.files <- lapply(class.names.files, function(x) fread(x, sep = "\t", data.table = FALSE))

## ------ removing mono-exonic transcripts ------
# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
class.files <- lapply(class.files, function(x) x %>% filter(!structural_category_exons %in% c("Intergenic_1","intergenic_1")))

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
mono.multi.class.files <- lapply(class.files, function(x) x %>% filter(monomulti == TRUE))
class.files <- lapply(class.files, function(x) x %>% filter(monomulti == FALSE))


## ------ filtering by expression ------

# filter isoforms that are only detected in the whole dataset, 2 reads, 2 samples
class.files$glob_SQ <- class.files$glob_targ_SQ %>% filter(whole_nsamples >= 2, whole_nreads >= 2)
class.files$targ_SQ <- class.files$glob_targ_SQ %>% filter(targeted_nsamples >= 2, targeted_nreads >= 2)

# known genes in whole dataset
class.files$glob_SQ_annoGene <- class.files$glob_SQ %>% filter(!grepl("novelGene", associated_gene)) %>% mutate(novelTranscript = ifelse(associated_transcript == "novel","Novel","Known"))

# annotated genic features
annoGenesStats <- list(
  novelTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel",],
  annoTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript != "novel",],
  NIC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NIC",],
  NNC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NNC",]
)

# subset by number of FL reads
femaleReads <- class.files$glob_SQ %>% select(all_of(femaleWhole)) %>% apply(., 1, sum)
maleReads <- class.files$glob_SQ %>% select(all_of(maleWhole)) %>% apply(., 1, sum)
postReads <- class.files$glob_SQ %>% select(all_of(postWhole)) %>% apply(., 1, sum)
preReads <- class.files$glob_SQ %>% select(all_of(preWhole)) %>% apply(., 1, sum)
class.files$glob_SQ <- class.files$glob_SQ %>% mutate(FReads = femaleReads, MReads = maleReads, preReads = preReads, postReads = postReads)
class.files$glob_SQ$DevStatus <- apply(class.files$glob_SQ, 1, function(x) identify_dataset_by_counts(x[["postReads"]], x[["preReads"]], "postnatal","prenatal"))

femaleMedianReads <- class.files$glob_SQ %>% select(all_of(femaleWhole)) %>% apply(., 1, median)
maleMedianReads <- class.files$glob_SQ %>% select(all_of(maleWhole)) %>% apply(., 1, median)
postMedianReads <- class.files$glob_SQ %>% select(all_of(postWhole)) %>% apply(., 1, median)
preMedianReads <- class.files$glob_SQ %>% select(all_of(preWhole)) %>% apply(., 1, median)
class.files$glob_SQ <- class.files$glob_SQ %>% mutate(FReads_median = femaleMedianReads, MReads_median = maleMedianReads, 
                                                      preReads_median = preMedianReads, postReads_median = postMedianReads)
save(class.files, file = paste0(dirnames$wholetarg_SQ,"sqantifiltered_monoexonicfiltered_2reads2samples.RData"))

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
  whole_sex = vroom(paste0(dirnames$DTE,"DESeq2_whole_sex_normSig.csv"),delim = ","),
  ns1 = read.csv(paste0(dirnames$DTE,"ONTX_6127_19264_whole_normAll.csv"), header = F),
  ns2 = read.csv(paste0(dirnames$DTE,"ONTX_6127_19264_whole_normAll.csv"), header = F)
)
colnames(Exp$ns1) <- colnames(Exp$whole_group)
colnames(Exp$ns2) <- colnames(Exp$whole_group)
Exp$whole_group <- rbind(Exp$whole_group,Exp$ns1,Exp$ns2)
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