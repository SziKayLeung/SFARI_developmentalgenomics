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
c(targeted_lengths1,targeted_lengths2,targeted_lengths3) <- plot_lengths(tpostnatal, tprenatal, tpfpostnatal, tpfprenatal)

targeted_lengths1/targeted_lengths2/targeted_lengths3+targeted_lengthslot_layout(guides = 'collect')

(targeted_lengths2+targeted_lengths1)/(targeted_lengths4+targeted_lengths3)

## whole
c(whole_lengths1,whole_lengths2,whole_lengths3) <- plot_lengths(wpostnatal, wprenatal, wpfpostnatal, wpfprenatal)

whole_lengths1/whole_lengths2/whole_lengths3+whole_lengthslot_layout(guides = 'collect')

(whole_lengths2+whole_lengths1)/(whole_lengths4+whole_lengths3)

# age distribution

c(age1, age2, age3, age4) <- ages(phenotype)

(age2+age1)/(age4+age3)
