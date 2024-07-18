## ---------- Script -----------------
##
## Purpose: perform differential isoform usage analysis on mouse rTg4510 ONT and Iso-Seq targeted datasets
## Adapted tappAS DIU analysis scripts
## https://github.com/ConesaLab/tappAS/blob/master/scripts/DIU.R
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- source functions -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_usage.R"))

library("data.table")
library("stringr")
library("dplyr")
library("ggplot2")

## ---------- input -----------------

dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/",
  root_sfari = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/",
  wholetarg_SQ = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/",
  testing = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/"
)

# classification files
#class.names.files = paste0(dirnames$testing, "SEPTIN4_reads.txt")
class.names.files = "/lustre/home/sl693/testingDIUCommand/SEPTIN4_classification.txt"
#class.files <- SQANTI_class_preparation(class.names.files,"nstandard")
class.files <- read.table(class.names.files,sep="\t",as.is=T)
row.names(class.files) <- class.files$isoform

# phenotype
phenotype <- list(
  WholeTargeted = read.csv(paste0(dirnames$root_sfari, "/12_deseq2/WholeTargetedphenotype.csv")) %>% mutate(time = age)
)
phenotype$WholeTargeted <- phenotype$WholeTargeted %>% mutate(col = paste0(sample,"_",group))

# DESEQ results
#Exp <- list(targeted = fread(paste0(dirnames$testing,"SEPTIN4_norm_targeted_sex.csv")))
Exp <- list(targeted = fread("/lustre/home/sl693/testingDIUCommand/SEPTIN4_normalised_expression.txt"))
Exp$targeted <- merge(Exp$targeted, phenotype$WholeTargeted, by="sample")
Exp$targeted <- merge(Exp$targeted, class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], by = "isoform", all.x = T)

# expression
ExpMod <- list(
  ontTargTranRaw = class.files %>% dplyr::select(associated_gene, contains("Targeted")) %>% filter(associated_gene == "SEPTIN4"),
  ontTargTranNorm = Exp$targeted %>%
    dplyr::select(sample,normalised_counts, isoform) %>% tidyr::spread(., sample, value = normalised_counts) %>%
    tibble::remove_rownames(.) %>% tibble::column_to_rownames(var="isoform")
)

# factors
Exp$targeted  <- list(
  targ_group = phenotype$WholeTargeted %>% filter((grepl("Targeted",sample))) %>% dplyr::select(sample, group) %>% magrittr::set_rownames(.$sample) %>%
    dplyr::rename("Replicate" = "group") %>%
    mutate(Replicate = ifelse(Replicate == "Prenatal",1,2))
)
factorsInput <- lapply(factorsInput, function(x) x[order(x$Replicate),, drop = FALSE])

resultsDIU <- list(
  ontTargGroup = runDIU(transMatrixRaw=ExpMod$ontTargTranRaw,transMatrix=ExpMod$ontTargTranNorm,classf=class.files,
                        myfactors=factorsInput$targ_group,filteringType="FOLD",filterFC=2)
)

IFSEPTIN4 <- plotIF(gene="SEPTIN4",ExpInput=ExpMod$ontTargTranRaw,pheno=phenotype$WholeTargeted,cfiles=class.files,design="case_control",majorIso=NULL,rank=4,isoSpecific=NULL)

plot_grid(plotlist=IFSEPTIN4)

Rscript ${LOGEN_ROOT}/differential_analysis/run_DIU_commandLine.R -f PRDM16-DT_classification.txt  -n PRDM16-DT_normalised_expression.txt

[sl693@mrc-comp085 2_differential_analysis]$ cat 3_targeted_dataset_output.R
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