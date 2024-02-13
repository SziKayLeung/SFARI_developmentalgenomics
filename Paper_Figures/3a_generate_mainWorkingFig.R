plot_grid(          plot_grid(fig3Diff$topRankedSex1Vis, fig3Diff$topRankedSex1,rel_heights = c(0.2,0.8),ncol=1),
          plot_grid(fig3Diff$topRankedSex2Vis, fig3Diff$topRankedSex2,rel_heights = c(0.2,0.8),ncol=1),
          plot_grid(fig3Diff$topRankedSex3Vis, fig3Diff$topRankedSex3,rel_heights = c(0.2,0.8),ncol=1),
          plot_grid(fig3Diff$topRankedSex4Vis, fig3Diff$topRankedSex4,rel_heights = c(0.2,0.8),ncol=1), nrow = 1)

plot_grid(
  plot_grid(Fig4Targeted$GRIA3DIU,Fig4Targeted$GRIA3DGE, nrow = 1),
  plot_grid(Fig4Targeted$GRIA3DTE1, Fig4Targeted$GRIA3DTE2, nrow = 1), rel_heights = c(0.6,0.4), nrow = 2)

function(gene, pathDIU){
  
  #Exp <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/targeted/group/",gene,"_normalised_expression.txt"),header=T)
  Exp <- read.csv(paste0(pathDIU,"/",gene,"_normalised_expression.txt"),header=T)
  Exp <- Exp %>% tidyr::spread(sample,normalised_counts) %>% tibble::column_to_rownames(var = "isoform")
  p <- plotIF(gene=gene,ExpInput=Exp,pheno=phenotype$WholeTargeted,cfiles=class.files$glob_targ_SQ,design="case_control",rank=5,majorIso=NULL)[[2]]
  return(p)
  
}


plotIFTargetedbyGene("CLIP3", paste0(dirnames$DIU,"whole/allSex"))


