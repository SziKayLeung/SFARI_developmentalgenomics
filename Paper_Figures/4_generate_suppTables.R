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

[sl693@mrc-comp083 Paper_Figures]$ cat 4b_generate_suppTables.R

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
plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$whole_group,gene="MECP2")
plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$targeted_group,gene="MECP2")

plot_grid(
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$whole_group,"group","MECP2") + ylim(2,3),
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$targeted_group,"group","MECP2") + ylim(2,3)
)

plot_grid(
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$whole_group,"group","EPN2"),
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$targeted_group,"group","EPN2")
)
[sl693@mrc-comp083 Paper_Figures]$ cat 5_ficle.R
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))


ficleDir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp"

pGRIK2 <- list(
  overview = plot_dendro_Tgene(ficleDir, "GRIK2"),
  A5A3 = plot_A5A3_Tgene(ficleDir, "GRIK2"),
  ES = plot_ES_Tgene(ficleDir,"GRIK2",class.files$glob_targ_SQ)[[1]],
  tracks = p
)

pTRIO <- list(
  overview = plot_dendro_Tgene(ficleDir, "TRIO"),
  A5A3 = plot_A5A3_Tgene(ficleDir, "TRIO"),
  ES = plot_ES_Tgene(ficleDir,"TRIO",class.files$glob_targ_SQ)[[1]],
  tracks = p2
)

SFARI_genes <- c("SHANK3","MECP2","CHD8","CNATNAP2","PTEN","SCN1A","SYNGAP1","ADNP","ARD1B","DYRK1A")
class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene %in% SFARI_genes,"isoform"]

WholeDESeqSig$sex %>% filter(!chrom %in% c("chrX","chrY")) %>%
  left_join(., class.files$glob_targ_SQ[,c("isoform","associated_gene")], by = "isoform") %>%
  mutate(SFARI_1_2 = ifelse(associated_gene %in% SFARI_CLASS_1_2,TRUE,FALSE),
         SFARI_1_2_S = ifelse(associated_gene %in% SFARI_CLASS_1_2_S ,TRUE,FALSE))

plot_grid(plot_trans_exp_individual("ONT5_6069_349",class.files$glob_targ_SQ,Exp$whole_sex,"sex") + facet_grid(~group),
          plot_trans_exp_lifetime("ONT5_6069_349",class.files$glob_targ_SQ,Exp$whole_sex)
)

plot_grid(
  plot_grid(pGRIK2[[1]],pGRIK2[[3]],pGRIK2[[2]] + theme(legend.position="top"),ncol=1),
  pGRIK2[[4]]
)

plot_grid(
  plot_grid(pTRIO[[1]],pTRIO[[3]] + scale_x_discrete(breaks = seq(10, 60, by = 10)),pTRIO[[2]] +
              scale_x_discrete(breaks = seq(10, 60, by = 10)) +
              theme(legend.position="top"),ncol=1),
  pTRIO[[4]], rel_widths = c(0.5,0.5)
)
[sl693@mrc-comp083 Paper_Figures]$ cat 5b_ficle.R
library(ggplot2)
library(ggtranscript)
library(plyranges)
library(vroom)

y<-class.files$glob_targ_SQ
x1 <- read_gff ("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/sqanti_with_ORF/WholeTargeted_with_ORF_corrected_renamed.gtf.cds.gff")
#x <- as.data.frame(x1 %>% filter(transcript_id %in% class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene == "CDH8","isoform"][1:10]))
dat2 <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp/TRIO/Stats/TRIO_Exonskipping_tab.csv")
#dat2 <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp/GRIK2/Stats/GRIK2_Exonskipping_tab.csv")
extensiveES <- dat2 %>% group_by(transcript_id) %>% tally() %>% filter(n > 40)
x <- as.data.frame(x1 %>% filter(transcript_id %in%
                                   c(as.character(extensiveES$transcript_id[1:30]),
                                     y[y$associated_gene == "TRIO" & y$structural_category == "FSM","isoform"][1:10],
                                     c("ONT5_1471_4273","ONT5_1471_4815","ONT5_1471_4942","ONT5_1471_1636"))))

#head(y[y$associated_gene == "TRIO",] %>% arrange(-length))

#x <- x[x$transcript_id %in% c(unique(x$transcript_id)[1:20],"ONT16_2235_233","ONT16_2235_493"),]
#cpat<-vroom('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.ORF_prob.best.tsv'))
#x<-read_gff('~/isoform_analysis/data/examples/ANKRD12.gff')
#x<-as.data.frame(x)
#x<-x[which(x$transcript_id%in%y$isoform),]
#x<-x[which(x$transcript_id%in%cpat[which(cpat$Coding_prob>0.364),]$seq_ID),]
#x<-x[which(x$transcript_id%in%tmp),]

p2 <- plot_CDS(x)
plot_CDS <- function(x){
  x_exons <- x %>% dplyr::filter(type == "exon")
  head(x_exons)
  
  # extract cds
  x_cds <- x %>% dplyr::filter(type == "CDS")
  head(x_cds)
  
  x_cds_w_stop <- x_cds %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      end = ifelse(end == max(end), end + 3, end)
    ) %>%
    dplyr::ungroup()
  
  # add_utr() adds ranges that represent the UTRs
  x_cds_utr <- add_utr(
    x_exons,
    x_cds_w_stop,
    group_var = "transcript_id"
  )
  x_cds_utr_rescaled <-
    shorten_gaps(
      exons = x_cds_utr,
      introns = to_intron(x_cds_utr, "transcript_id"),
      group_var = "transcript_id"
    )
  
  z<-read_gff('/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/0_output/gencode.v44.annotation.TRIO.gff3')
  z<-as.data.frame(z)
  z<-z[which(z$gene_name=='TRIO'),]
  
  z_exons <- z %>% dplyr::filter(type == "exon")
  z_exons %>% head()
  
  # extract cds
  z_cds <- z %>% dplyr::filter(type == "CDS")
  z_cds %>% head()
  
  z_cds_w_stop <- z_cds %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      end = ifelse(end == max(end), end + 3, end)
    ) %>%
    dplyr::ungroup()
  
  # add_utr() adds ranges that represent the UTRs
  z_cds_utr <- add_utr(
    z_exons,
    z_cds_w_stop,
    group_var = "transcript_id"
  )
  z_cds_utr_rescaled <-
    shorten_gaps(
      exons = z_cds_utr,
      introns = to_intron(z_cds_utr, "transcript_id"),
      group_var = "transcript_id"
    )
  
  
  x<-rbind(x, z[colnames(x)])
  
  p <- x_cds_utr_rescaled %>%
    dplyr::filter(type == "CDS") %>%
    ggplot(., aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(fill='blue') +
    geom_range(
      data = x_cds_utr_rescaled %>% dplyr::filter(type == "UTR"),
      height = 0.25,
      fill = "white"
    ) +
    geom_intron(
      data = to_intron(
        x_cds_utr_rescaled %>% dplyr::filter(type != "intron"),
        "transcript_id"
      ),
      arrow.min.intron.length = 110
    ) +
    geom_range(data=z_cds_utr_rescaled %>% dplyr::filter(type == "CDS"),
               aes(
                 xstart = start,
                 xend = end,
                 y = transcript_id
               ),fill='blue') +
    geom_range(
      data = z_cds_utr_rescaled %>% dplyr::filter(type == "UTR"),
      height = 0.25,
      fill = "white"
    ) +
    geom_intron(
      data = to_intron(
        z_cds_utr_rescaled %>% dplyr::filter(type != "intron"),
        "transcript_id"
      ),
      arrow.min.intron.length = 110
    ) + theme_classic() + labs(x=NULL,y="Transcripts")
  
  
  return(p)
}