#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary tables for isoform developmental paper
##
## ---------------------------------
## fig1:
## fig2: 
## fig3:
## fig4: Targeted dataset


library(xlsx)
SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs/SuppTables")


##************** Whole DESeq 
# DTE 
DESeqCols <- c("isoform","chrom", "associated_gene", "associated_transcript","log2FoldChange","lfcSE","pvalue","padj","structural_category","subcategory")
DESeqColsRenamed <- c("Isoform","Chromosome", "Associated gene", "Associated transcript","Log2FC","lfcSE", "p-value","FDR","Structural category","Subcategory")
WholeDESeqSigOut <- lapply(WholeDESeqSig, function(x) x %>% arrange(padj) %>% select(all_of(DESeqCols)) %>% `colnames<-`(DESeqColsRenamed))
write.csv(WholeDESeqSigOut$sex,paste0(output_dir,"/WholeDTESex.csv"))
write.csv(WholeDESeqSigOut$age,paste0(output_dir,"/WholeDTEGroup.csv"))
# DGE 
DESeqGeneCols <- c("associated_gene","log2FoldChange","lfcSE","pvalue","padj")
DESeqGeneColsRenamed <- c("Associated gene","Log2FC","lfcSE", "p-value","FDR")
WholeDESeqGeneSigOut <- lapply(WholeDESeqGeneSig, function(x) x %>% mutate(TargetGene = ifelse(associated_gene %in% TargetGene,TRUE,FALSE)))
WholeDESeqGeneSigOut <- lapply(WholeDESeqGeneSigOut, function(x) x %>% arrange(padj) %>% select(all_of(c(DESeqGeneCols,"TargetGene"))) %>% `colnames<-`(c(DESeqGeneColsRenamed,"TargetGene")))
write.csv(WholeDESeqGeneSigOut$sex,paste0(output_dir,"/WholeDGESex.csv"))
write.csv(WholeDESeqGeneSigOut$age,paste0(output_dir,"/WholeDGEGroup.csv"))

# Targeted DESeq
# DTE (age only, no sig for sex)
TargetedDESeq2SigOut <- TargetedDESeqSig$age %>% arrange(padj) %>% select(all_of(DESeqCols[ !DESeqCols == 'chrom'])) %>% `colnames<-`(DESeqColsRenamed[ !DESeqColsRenamed == 'Chromosome'])
write.csv(TargetedDESeq2SigOut,paste0(output_dir,"/TargetedDTEGroup.csv"))
# DGE 
TargetedDESeqGeneSigOut <- lapply(TargetedDESeqGeneSig, function(x) x %>% arrange(padj) %>% select(all_of(DESeqGeneCols)) %>% `colnames<-`(DESeqGeneColsRenamed))
write.csv(TargetedDESeqGeneSigOut$sex,paste0(output_dir,"/TargetedDGESex.csv"))
write.csv(TargetedDESeqGeneSigOut$age,paste0(output_dir,"/TargetedDGEGroup.csv"))
write.table(WholeDESeq2Sig$sex,paste0(dirnames$output,"anno_whole_sex_nomonointergenic_Sig.csv"),quote = F,row.names = F,col.names = T)
write.table(WholeDESeq2Sig$age,paste0(dirnames$output,"anno_whole_group_nomonointergenic_Sig.csv"),quote = F,row.names = F,col.names = T)
# DIU
DIUColsRenamed <- c("Associated gene","p-value","FDR", "Podium change","Total change")
DIUSigOut <- lapply(DIUSig, function(x) x %>% arrange(FDR) %>% `colnames<-`(DIUColsRenamed))
write.csv(DIUSigOut$targetedSex,paste0(output_dir,"/TargetedDIUSex.csv"))
write.csv(DIUSigOut$targetedAge,paste0(output_dir,"/TargetedDIUGroup.csv"))

