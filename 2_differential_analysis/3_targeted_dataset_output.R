library("dplyr")
DIU_targeted$FDR <- as.numeric(as.character(DIU_targeted$FDR))
DIU_targeted$p.value <- as.numeric(as.character(DIU_targeted$p.value))
DIUSig <- DIU_targeted %>% filter(FDR < 0.05) %>% arrange(FDR)
nrow(DIU_targeted)
nrow(DIUSig)
nrow(DIUSig %>% filter(podiumChange == "TRUE"))

LOGEN_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "differential_analysis/plot_usage.R"))

gene = "TRPC4"

phenotype$WholeTargeted <- phenotype$WholeTargeted %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), col = paste0(sample,"_",group))

plotIFbyGene <- function(gene){

  Exp <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/targeted/",gene,"_normalised_expression.txt"),header=T)
  Exp <- Exp %>% tidyr::spread(sample,normalised_counts) %>% tibble::column_to_rownames(var = "isoform")
  p <- plotIF(gene=gene,ExpInput=Exp,pheno=phenotype$WholeTargeted,cfiles=class.files$glob_targ_SQ,design="case_control",rank=5,majorIso=NULL)[[2]]
  return(p)

}

DIUPlots <- list()
for(i in DIUSig$Gene){
  DIUPlots[[i]] <- plotIFbyGene(i)
}

pdf("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/Paper_Figures/outputFigs/DIUTargeted.pdf")
DIUPlots
dev.off()

length(intersect(TargetedDESeqSig$age$associated_gene,DIUSig$Gene))
TargetedDESeqSigNumTranscripts <- TargetedDESeqSig$age %>% group_by(associated_gene) %>% tally()
TargetedDESeqSigNumTranscripts %>% filter(associated_gene %in% DIUSig$Gene)
TargetedDESeqSig$age %>% filter(associated_gene == "CD38")

plot_trans_exp_individual("ONT11_2781_9724",class.files$glob_targ_SQ,Exp$targeted_group,"group")
plot_trans_exp_individual("ONT4_685_5082",class.files$glob_targ_SQ,Exp$targeted_group,"group")