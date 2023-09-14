## ---------- Script -----------------
##
## Purpose: perform differential analysis on pancreas targeted datasets using linear regression
## Transcript level separate analysis
## Gene level separate analysis
##
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
##
## --------- Notes -----------------
## Expression


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))


## ---------- source functions -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))


## ---------- input -----------------

# directory names
dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/",
  sqanti = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/1_SQANTI3Filtered/"
)

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$root, "phenotype.csv"), 
  expression = paste0(dirnames$root, "demux_fl_count.csv"),
  classfiles = paste0(dirnames$sqanti, "merged_RulesFilter_result_classification.targetgenes_counts.txt")
)


input <- list()
input$ontPhenotype <- read.csv(input_files$ontPhenotype, header = T) 
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$expression <- read.csv(input_files$expression) %>% .[.$isoform != "0",] %>% filter(isoform %in% input$classfiles$isoform)
colnames(input$ontPhenotype) <- c("sample","time","group","age")

## ---------- ONT: Creating DESeq2 object and analysis -----------------

# run DESeq2
ontResTran <- list(
  wald = run_DESeq2(test="Wald",input$expression,input$ontPhenotype,threshold=10,exprowname="isoform",controlname="M",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$expression,input$ontPhenotype,threshold=10,exprowname="isoform",controlname="M",design="time_series",interaction="On"),
)

ontResTranAnno <- lapply(ontResTran, function(x) anno_DESeq2(x,input$classfiles,input$ontPhenotype,controlname="CONTROL",level="transcript",sig=0.1))

# split by genotype and age effects
ontResTranEffects <- do.call(rbind,dissect_DESeq2(wald=ontResTranAnno$wald$anno_res,lrt=ontResTranAnno$lrt$anno_res))
ontResTranEffects

## ---------- Iso-Seq: Creating DESeq2 object and analysis -----------------

# run DESeq2
isoResTran <- list(
  wald = run_DESeq2(test="Wald",input$isoExpression,input$isoPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$isoExpression,input$isoPhenotype,threshold=10,exprowname="isoform",controlname="CONTROL",design="time_series",interaction="On"),
  wald8mos = run_DESeq2(input$isoExpression,input$isomos8Phenotype,exprowname="isoform",controlname="CONTROL",design="case_control",interaction="On",test="Wald")
)

isoResTranAnno <- lapply(isoResTran, function(x) anno_DESeq2(x,input$classfiles,input$isoPhenotype,controlname="CONTROL",level="transcript",sig=0.1))

# split by genotype and age effects
isoResTranEffects <- do.call(rbind,dissect_DESeq2(wald=isoResTranAnno$wald$anno_res,lrt=isoResTranAnno$lrt$anno_res))
isoResTranEffects


## ---------- Iso-Seq Differential gene expression -----------------

# run DESeq2
ontResGene <- list(
  wald = run_DESeq2(test="Wald",input$gene_expression %>% select(-associated_gene),input$ontPhenotype,threshold=10,controlname="CONTROL",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$gene_expression %>% select(-associated_gene),input$ontPhenotype,threshold=10,controlname="CONTROL",design="time_series",interaction="On")
)

isoResGene <- list(
  wald = run_DESeq2(test="Wald",input$gene_expression %>% select(-associated_gene),input$isoPhenotype,threshold=10,controlname="CONTROL",design="time_series",interaction="On"),
  lrt = run_DESeq2(test="LRT",input$gene_expression %>% select(-associated_gene),input$isoPhenotype,threshold=10,controlname="CONTROL",design="time_series",interaction="On")
)

# annotate results
ontResGeneAnno <- lapply(ontResGene, function(x) anno_DESeq2(x,input$classfiles,input$ontPhenotype,controlname="CONTROL",level="gene",sig=0.1))
isoResGeneAnno <- lapply(isoResGene, function(x) anno_DESeq2(x,input$classfiles,input$isoPhenotype,controlname="CONTROL",level="gene",sig=0.1))


## ---------- Output -----------------

saveRDS(ontResTranAnno, file = paste0(dirnames$output, "/Ont_DESeq2TranscriptLevel.RDS"))
saveRDS(isoResTranAnno, file = paste0(dirnames$output, "/IsoSeq_DESeq2TranscriptLevel.RDS"))
saveRDS(ontResGeneAnno, file = paste0(dirnames$output, "/Ont_DESeq2GeneLevel.RDS"))
saveRDS(isoResGeneAnno, file = paste0(dirnames$output, "/IsoSeq_DESeq2GeneLevel.RDS"))