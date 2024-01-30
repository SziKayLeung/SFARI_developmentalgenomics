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
  
  output = paste0(root_sfari,"/0_output/"),
  
  diu = paste0(root_dir, "RBFetal/5_diu/")
)

TargetGene = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/vc362/fetalBrain/genes.txt")[["V1"]]
# list of SZ and ASD genes in targeted panel
TargetGeneSZASD = read.csv(paste0(root_sfari,"0_metadata/TargetGeneByDisease.csv"))

## -------------- Final classification files -------------

# class.files
# list: glob_targ_SQ, glob_SQ, targ_SQ
load(file = paste0(dirnames$wholetarg_SQ,"all_filtered_classification_2reads2samples_noMonoIntergenic.RData"))

wholesamples <- colnames(class.files$glob_targ_SQ)[grepl("Whole", colnames(class.files$glob_targ_SQ))]
targetedsamples <- colnames(class.files$glob_targ_SQ)[grepl("Targeted", colnames(class.files$glob_targ_SQ))]
matchedsamples <- intersect(gsub("^.*?Whole","",wholesamples),gsub("^.*?Targeted","",targetedsamples))
targetedmatchedsamples <- paste0("Targeted",matchedsamples)
wholematchedsamples <- paste0("Whole",matchedsamples)


## -------------- DESeq2 ----------------

# TargetedDESeq (all results)
load(file = paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_removinglateprenatal.RData"))
# TargetedDESeqSig (sig FDR < 0.05)
load(file = paste0(root_dir,"RBFetal/4_deseq2/anno_targeted_removinglateprenatal_sig.RData"))

# WholeDESeq
load(file = paste0(root_dir,"RBFetal/4_deseq2/anno_whole_removinglateprenatal_TEST.RData"))
# WholeDESeqSig (sig FDR < 0.05)
load(file = paste0(root_dir,"RBFetal/4_deseq2/anno_whole_removinglateprenatal_TEST_sig.RData"))


## ------------- Phenotype files -------------------
phenotype <- list(
  WholeTargeted = fread(paste0(root_dir, "RBFetal/00_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
)
phenotype$WholeTargeted <- phenotype$WholeTargeted %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), col = paste0(sample,"_",group))


## -------------- differenetial gene expression -------------
WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_prenatalsex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_whole_gene_group_removinglateprenatal.csv")))
)

TargetedDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_gene_prenatalsex_removinglateprenatal.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"RBFetal/4_deseq2/filtered_output_targeted_gene_group_removinglateprenatal.csv")))
)


## -------------- differential isoform usage ----------------

read_DIU <- function(inputPath){
  DIU_targeted <- list.files(path=inputPath,full.names = T, pattern = "resultDIU")
  if(length(DIU_targeted) > 0){
    DIU_targeted <- lapply(DIU_targeted,function(x) read.table(x)[-1,])
    DIU_targeted <- do.call(rbind, DIU_targeted)
    colnames(DIU_targeted) <- c("Gene","p.value","FDR","podiumChange","totalChange")
    DIU_targeted <- DIU_targeted %>% mutate(FDR = as.numeric(as.character(FDR)))
  }else{
    return(NULL)
  }
}

DIU <- list(
  targetedAge = read_DIU(paste0(dirnames$diu,"targeted/group")),
  targetedSex = read_DIU(paste0(dirnames$diu,"targeted/sex")),
  wholeAge = read_DIU(paste0(dirnames$diu,"whole/group")),
  wholeSex = read_DIU(paste0(dirnames$diu,"whole/sex"))
)
DIUSig <- lapply(DIU, function(x) x[x$FDR <= 0.05, ])


## -------------- normalized counts ----------------

#Exp
load(file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_removinglateprenatal.RData"))
load(file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_gene_removinglateprenatal.RData"))
load(file = paste0(root_dir,"RBFetal/4_deseq2/filtered_output_norm_gene_sig_removinglateprenatal.RData"))
Exp <- Exp2

## -------------- cpat ----------------

Cpat <- list(
  whole = data.table::fread(paste0(root_dir,"RBFetal/2_cpat_tc20bp/Whole.ORF_prob.best.tsv"))
)
#save(Cpat, file = paste0(root_dir,"RBFetal/2_cpat_tc20bp/Whole.ORF_prob.best.RData"))
# keep  only the list of isoforms in the final dataset
Cpat$whole <- Cpat$whole %>% filter(seq_ID %in% class.files$glob_SQ$isoform)

## -------------- gtf ----------------

gtf <- list(
  glob_targ = rtracklayer::import(paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf"))
  #ref = rtracklayer::import(paste0(dirnames$output,"gencode.v40.ggTransRefGenes.gtf"))
  #ref = rtracklayer::import(paste0(dirnames$output,"refExons.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))
gtf$ref <- data.table::fread(paste0(dirnames$output,"refExons.gtf")) %>% dplyr::rename("gene_id" = "gene_name") %>% mutate(type = "exon")
gtf$merged <- rbind(gtf$glob_targ[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                    gtf$ref[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

GI <- c("GRIN2A","GRIA3","SEPTIN4","RTN4","MBP","RPS4Y1","XIST","ADD3","CNTNAP2","ANKRD12")
RefIsoforms <- lapply(GI, function(x) unique(gtf$ref[gtf$ref$gene_id == x & !is.na(gtf$ref$transcript_id), "transcript_id"]))
names(RefIsoforms ) <-GI


## -------------- disease list ----------------

diseasegenelists <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/0_utilities/disease_list/"
SFARI <- read.csv(paste0(diseasegenelists,"SFARI-Gene_genes_07-17-2023release_09-26-2023export.csv"),header=T)
SFARI_CLASS_1_2_S <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2 |SFARI$syndromic == 1)
SFARI_CLASS_1_2 <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2)
disease_list <- list(
  SCHEMA = read.table(paste0(diseasegenelists,"SCHEMA_Oct2023.csv"), sep=",", header=T, stringsAsFactors = F),
  DDG2P = read.table(paste0(diseasegenelists,"DDG2P_26_9_2023.csv"), sep=",", header=T, stringsAsFactors = F) %>% filter(confidence.category != "limited")
)
