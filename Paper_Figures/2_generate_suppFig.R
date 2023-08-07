#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary figures
##
## ---------------------------------


# sensitivity curve of filtering in targeted dataset
no_of_isoforms_sample(class.files$targ_SQ)

plot_cupcake_collapse_sensitivity(class.files$targ_SQ,"All target genes")
