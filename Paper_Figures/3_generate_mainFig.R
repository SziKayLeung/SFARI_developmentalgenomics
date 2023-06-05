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


fig4Targeted <- list(
  # number of isoforms 
  numIso = numIso(class.files$targ_SQ_fil),
  # number of isoforms by category
  numIsoCate = numIsoCate(class.files$targ_SQ_fil)
)