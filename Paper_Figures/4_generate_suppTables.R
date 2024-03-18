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

WholeDTENumGene <- WholeDTE$age %>% group_by(associated_gene, dirAcrossDev) %>% tally() %>% 
  as.data.frame() %>% reshape(., idvar = "associated_gene", timevar = "dirAcrossDev", direction = "wide")
WholeDTENumGene[is.na(WholeDTENumGene)] <- 0
write.csv(WholeDTENumGene,paste0(dirnames$output,"/WholeDTENumGene.csv"), quote = F)

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

# number of transcripts pre- and post-natal
class.files$glob_SQ_annoGene %>% group_by(structural_category,DevStatus) %>% tally(name = "num") %>% 
  ggplot(., aes(x = structural_category, y = num, fill = DevStatus)) + geom_bar(stat = "identity", position="dodge") +
  labs(x = "Structural category", y = "Number of transcripts of annotated genes") +
  scale_fill_discrete(label = c("Both","Post-natal only", "Pre-natal only"), name = NULL) +
  theme_classic() + theme(legend.position = "top")

uniquePrenatalTranscripts <- class.files$glob_SQ_annoGene %>% filter(DevStatus == "prenatal")
uniquePrenatalGenes <- unique(class.files$glob_SQ_annoGene %>% filter(DevStatus == "prenatal") %>% select(associated_gene))
write.table(uniquePrenatalGenes, paste0(dirnames$output,"uniquePrenatalGenes.txt"), row.names = F, col.names = F, quote = F)

# ratio of prenatal vs postnatal
TallyNumNovelKnown <- as.data.frame(class.files$glob_SQ_annoGene %>% dplyr::group_by(novelTranscript, associated_gene) %>% tally())
TallyNumNovelKnown <- reshape(TallyNumNovelKnown, idvar = "associated_gene", timevar = "novelTranscript", direction = "wide")
geneNum <- Reduce(function(...) merge(..., all=T, by = "associated_gene"), 
                  list(class.files$glob_SQ_annoGene %>% dplyr::group_by(associated_gene) %>% tally(name = "totalNumTranscripts"),
                       TallyNumNovelKnown,
                       class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,] %>% group_by(associated_gene) %>% tally(name = "prenatalTranscripts"),
                       class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,] %>% group_by(associated_gene) %>% tally(name = "postnatalTranscripts")))
geneNum  <- geneNum %>% mutate(ratioPrevsPost = prenatalTranscripts/postnatalTranscripts, ratioPostvsPre = postnatalTranscripts/prenatalTranscripts)
geneNum  <- geneNum %>% mutate(DiffGeneExp = ifelse(associated_gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE))
WholeDESeqGeneSig$age$associated_gene
colnames(geneNum)[1:4] <- c("associated_gene","totalTranscripts","totalNovelTranscripts","totalKnownTranscripts")
write.csv(geneNum, paste0(dirnames$output,"NumberofTranscriptsPrevsPost.csv"), quote=F, row.names = F)

write.csv(DIUSig$wholeAllAge, paste0(dirnames$output,"DIUSigWholeDevelopment.csv"), quote=F, row.names = F)
write.csv(DIUSig$wholeAllSex, paste0(dirnames$output,"DIUSigWholeSex.csv"), quote=F, row.names = F)

# SQANTI output table
write.table(class.files$glob_SQ, paste0(dirnames$output,"SQANTI_WholeTranscriptomeDataset_Finalized.txt"), quote = F, sep  = "\t")

# description of target genes
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs/SuppTables")
dat <- data.frame()

for(i in 1:length(TargetGene)){
  gene = as.character(TargetGene[i])
  dat[i,1] <- gene
  dat[i,2] <- ifelse(nrow(WholeDESeqGeneSig$age[WholeDESeqGeneSig$age$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,3] <- ifelse(nrow(TargetedDESeqGeneSig$age[TargetedDESeqGeneSig$age$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,4] <- ifelse(nrow(WholeDESeqGeneSig$sex[WholeDESeqGeneSig$sex$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,5] <- ifelse(nrow(TargetedDESeqGeneSig$sex[TargetedDESeqGeneSig$sex$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,6] <- nrow(WholeDESeqSig$age %>% filter(associated_gene == gene))
  dat[i,7] <- nrow(TargetedDESeqSig$age %>% filter(associated_gene == gene))
  dat[i,8] <- nrow(WholeDESeqSig$sex %>% filter(associated_gene == gene))
  dat[i,9] <- nrow(TargetedDESeqSig$sex %>% filter(associated_gene == gene))
}
colnames(dat) <- c("TargetGene","WholeDGEGroup","TargetedDGEGroup","WholeDGESex","TargetedDGESex","NumWholeDTEGroup","NumTargetedDTEGroup","NumWholeDTESex","NumTargetedDTESex")

dat <- dat %>% mutate(DGEGroup = ifelse(WholeDGEGroup == TargetedDGEGroup, TRUE,FALSE),
                      DGESex = ifelse(WholeDGESex == TargetedDGESex, TRUE,FALSE))


write.csv(dat,paste0(output_dir,"/TargetGenesDEAnalysis.csv"))
