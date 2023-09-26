#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary figures
##
## ---------------------------------
library(patchwork)

# sensitivity curve of filtering in targeted dataset
no_of_isoforms_sample(class.files$targ_SQ)

plot_cupcake_collapse_sensitivity(class.files$targ_SQ,"All target genes")

# read lengths
## targeted
tpostnatal$V3<-'Pre-filter'
tprenatal$V3<-'Pre-filter'
tpfprenatal$V3<-'Post-filter'
tpfpostnatal$V3<-'Post-filter'

targeted_lengths <- plot_lengths(bind_rows(tpostnatal, tprenatal, tpfpostnatal, tpfprenatal))
targeted_lengths

## whole
wpostnatal$V3<-'Pre-filter'
wprenatal$V3<-'Pre-filter'
wpfprenatal$V3<-'Post-filter'
wpfpostnatal$V3<-'Post-filter'
whole_lengths <- plot_lengths(bind_rows(wpostnatal, wprenatal, wpfpostnatal, wpfprenatal))
whole_lengths

# age distribution

c(age1, age2, age3, age4) <- ages(phenotype)

(age2+age1)/(age4+age3)
