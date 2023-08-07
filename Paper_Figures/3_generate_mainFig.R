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
output_dir = paste0(SC_ROOT,"/figures/")

fig1Whole <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$glob_SQ, input_title = "", glimit = 10),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$glob_SQ)
)


fig2Targeted <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$targ_SQ_fil, input_title = "", glimit = 10),
  # number of genes
  numIsogene = numIsoGene(class.files$targ_SQ_fil),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$targ_SQ_fil)
  
)

pdf(paste0(output_dir,"/TargetedFigures.pdf"), width = 18, height = 12)
plot_grid(
  plot_grid(fig2Targeted$numIsogene, fig2Targeted$numIso, labels = c("D","E"), rel_widths = c(0.4,0.6),label_size = 20),
  plot_grid(NULL, fig2Targeted$numIsoCate, labels = c("F"),rel_widths = c(0.3,0.7),label_size = 20),
  nrow=2
)
dev.off()