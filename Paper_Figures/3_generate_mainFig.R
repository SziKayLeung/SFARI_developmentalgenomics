#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for main figures for isoform developmental paper
##
## ---------------------------------
## fig1:
## fig2: 
## fig3:
## fig4: Targeted dataset


SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs")


## ------ Figure 1: Whole dataset descriptives ------ 

fig1Whole <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$glob_SQ, input_title = "", glimit = 10),
  # number of genes
  numIsogene = numIsoGene(class.files$glob_SQ),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$glob_SQ)
)


## ------ Figure 2: Targeted dataset descriptives ------ 

fig2Targeted <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$targ_SQ, input_title = "", glimit = 10),
  # number of genes
  numIsogene = numIsoGene(class.files$targ_SQ),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$targ_SQ)
)


## ------ Figure 3: Differential expression analysis ------ 

plot_grid(fig1Whole$numIso,fig2Targeted$numIso, labels = c("A","B"))
plot_grid(fig1Whole$numIsogene,fig2Targeted$numIsogene, labels = c("A","B"))
plot_grid(fig1Whole$numIsoCate,fig2Targeted$numIsoCate, labels = c("A","B"))

# targeted dataset
Top10TargetedDESeq2Sex = plot_top_results(TargetedDESeq2$sex,Exp$targeted,"sex")
Top10TargetedDESeq2Age = plot_top_results(TargetedDESeq2$age,Exp$targeted,"group")

# whole transcriptome
Top10WholeDESeq2Sex = plot_top_results(WholeDESeq2$sex,Exp$wholeanno,"sex")
Top10WholeDESeq2Age = plot_top_results(WholeDESeq2$age,Exp$wholeanno,"group")


message("Number of significant differentially expressed transcripts across post and pre-natal: ", nrow(WholeDESeq2$age %>% filter(padj < 0.05)))
message("Number of significant differentially expressed transcripts across sex: ", nrow(WholeDESeq2$sex %>% filter(padj < 0.05)))


## ------ Figure 4: CPAT, alternative splicinge events ------ 

# tracks
ggTranPlots(inputgtf=gtf$glob_targ,classfiles=class.files$glob_targ_SQ,isoList=c("ONT1_2_13","ONT1_2_41"),colours = c("red","blue"),simple=TRUE)


## ------ Output ------ 

pdf(paste0(output_dir,"/MainFigure1.pdf"), width = 18, height = 12)
plot_grid(
  plot_grid(fig1Whole$numIsogene, fig1Whole$numIso, labels = c("D","E"), rel_widths = c(0.3,0.7),label_size = 20),
  plot_grid(NULL, fig1Whole$numIsoCate, labels = c("F"),rel_widths = c(0.3,0.7),label_size = 20),
  nrow=2
)
dev.off()

pdf(paste0(output_dir,"/MainFigure2.pdf"), width = 18, height = 12)
plot_grid(
  plot_grid(fig2Targeted$numIsogene, fig2Targeted$numIso, labels = c("D","E"), rel_widths = c(0.4,0.6),label_size = 20),
  plot_grid(NULL, fig2Targeted$numIsoCate, labels = c("F"),rel_widths = c(0.3,0.7),label_size = 20),
  nrow=2
)
dev.off()



plot_volcano(diff_results=WholeDESeq2$sex)
ggsave(file="WholeDeSeq2Sex.png", dpi=400)

plot_volcano(diff_results=WholeDESeq2$age)
ggsave(file="WholeDeSeq2Age.png", dpi=400)

plot_volcano(diff_results=TargetedDESeq2$sex)
ggsave(file="TargetedDeSeq2Sex.png", dpi=400)

plot_volcano(diff_results=TargetedDESeq2$age)
ggsave(file="TargetedDeSeq2Age.png", dpi=400)

