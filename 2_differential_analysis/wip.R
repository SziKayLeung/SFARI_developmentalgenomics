DIUSig$targetedAge
TargetedDESeqGeneSig$age
TargetedDESeqSig$age

TargetedAgeGenes <- list(unique(DIUSig$targetedAge$Gene),unique(TargetedDESeqGeneSig$age$associated_gene),unique(TargetedDESeqSig$age$associated_gene))
TargetedSexGenes <- list(unique(DIUSig$targetedSex$Gene),unique(TargetedDESeqGeneSig$sex$associated_gene),unique(TargetedDESeqSig$sex$associated_gene))

vennDEA <- function(listGenes,listnames){
  p <- plot_grid(venn.diagram(
    x = listGenes,
    category.names = listnames,
    filename = NULL,
    output = FALSE ,
    imagetype="png" ,
    compression = "lzw",
    lwd = 1,
    col=c("#440154ff", '#21908dff', '#fde725ff'),
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
    rotation = 1
  ))
  return(p)
}

vennDEA(TargetedAgeGenes,c("DIU","DGE","DTE"))
vennDEA(TargetedSexGenes,c("DIU","DGE","DTE"))
TargetedGenesDIUDTE <- setdiff(intersect(unique(DIUSig$targetedAge$Gene),unique(TargetedDESeqSig$age$associated_gene)),unique(TargetedDESeqGeneSig$age$associated_gene))


plotIFTargetedbyGene <- function(gene){
  
  Exp <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/targeted/group/",gene,"_normalised_expression.txt"),header=T)
  Exp <- Exp %>% tidyr::spread(sample,normalised_counts) %>% tibble::column_to_rownames(var = "isoform")
  p <- plotIF(gene=gene,ExpInput=Exp,pheno=phenotype$WholeTargeted,cfiles=class.files$glob_targ_SQ,design="case_control",rank=5,majorIso=NULL)[[2]]
  return(p)
  
}


plotDEAGene <- function(gene,transcript=NULL){
  print(gene)
  df <- TargetedDESeqSig$age[TargetedDESeqSig$age$associated_gene == gene,] %>% arrange(padj)
  if(is.null(transcript)){
    transcript = df$isoform
  }
  p <- list(
    DTE = plot_trans_exp_lifetime(transcript,class.files$glob_targ_SQ,Exp$targeted_group),
    DGE = plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$targeted_group,gene=gene),
    DIU = plotIFTargetedbyGene(gene)  
  )
  plist <- plot_grid(plotlist = p, nrow=1)
  return(list(plist,p))
}

OneDTEDIU <- TargetedDESeqSig$age %>% filter(associated_gene %in% TargetedGenesDIUDTE) %>% group_by(associated_gene) %>% tally() %>% .[.$n == 1, "associated_gene"]

pDEA <- lapply(OneDTEDIU$associated_gene, function(x) plotDEAGene(x))
names(pDEA) <- OneDTEDIU$associated_gene

GRIA3a <- plotDEAGene("GRIA3","ONTX_7115_9753")
GRIA3b <- plotDEAGene("GRIA3","ONTX_7115_973")

plot_grid(GRIA3a[[2]][[1]],GRIA3b[[2]][[1]],GRIA3b[[2]][[2]],GRIA3b[[2]][[3]])
nrow(TargetedDESeqSig$age[TargetedDESeqSig$age$associated_gene == "GRIA3",])


