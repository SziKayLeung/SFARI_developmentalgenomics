#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary figures
##
## ---------------------------------


SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs")


# sensitivity curve of filtering in targeted dataset
no_of_isoforms_sample(class.files$targ_SQ)

# comparison of whole vs targeted datasets across matched samples
comp = whole_vs_targeted_plots(class.files$glob_targ_SQ, wholematchedsamples, targetedmatchedsamples, TargetGene)

plot_grid(plotlist = comp[1:6], labels = c("A","B","C","D","E","F"))

plot_cupcake_collapse_sensitivity(class.files$targ_SQ,"All target genes")

# sensitivity plots
#pSensitivity(class.files$targ_SQ)
#ggsave(file=paste0(output_dir,"/cumulativeSensitivityTargeted.png"), dpi=400, width = 20, height = 20, units = "cm")

#pSensitivity(class.files$glob_SQ)
#ggsave(file=paste0(output_dir,"/cumulativeSensitivityWhole.png"), dpi=400, width = 20, height = 20, units = "cm")

pIF <- list(
  ontNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_ont$normAll,
                                                      pheno=phenotype$targeted_rTg4510_ont,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))),
  
  isoNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_iso$normAll,
                                                      pheno=phenotype$targeted_rTg4510_iso,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$isoDIUGeno$keptIso)))
)
for(i in 1:length(pIF)){names(pIF[[i]]) <- Targeted$Genes}

