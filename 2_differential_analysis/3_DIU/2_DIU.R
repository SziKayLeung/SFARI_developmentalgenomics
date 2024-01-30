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


## ---------- input -----------------

# directory names
dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/",
  root_sfari = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/",
  wholetarg_SQ = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/",
  meta = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/metadata"
)


TargetGene = read.table(paste0(dirnames$root_sfari, "/0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]
TargetGene20 = read.csv(paste0(dirnames$root_sfari, "/0_metadata/20SexDifferenceTargetGenes.csv"), header = F)[["V1"]]

# classification files
class.names.files <- list(
  
  # targeted + whole SQANTI dataset futher filtered by minimum 2 reads and 2 samples
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt"),
  
  # whole transcriptome SQANTI dataset further filtered by minimum 2 reads and 2 samples
  glob_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples.txt"),
  
  # targeted SQANTI dataset further filtered by minimum 2 reads and 2 samples
  targ_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples.txt")
  
)

class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# filter targeted dataset to just the target genes
class.files$targ_SQ <- class.files$targ_SQ %>% filter(associated_gene %in% TargetGene)


# DESEQ results
Exp <- list(
  targeted = fread(paste0(dirnames$root,"4_deseq2/output_norm_targeted_sex.csv"))
)
Exp$targeted <- merge(Exp$targeted, phenotype$WholeTargeted, by="sample")
Exp$targeted <- merge(Exp$targeted, class.files$glob_targ_SQ[,c("isoform","associated_gene","associated_transcript","structural_category")], by = "isoform", all.x = T)


# expression
ExpMod <- list(
  ontTargTranRaw = class.files$targ_SQ %>% dplyr::select(associated_gene, contains("Targeted")),
  ontTargTranNorm = Exp$targeted %>%
    dplyr::select(sample,associated_gene,normalised_counts, isoform) %>% spread(., sample, value = normalised_counts) %>%
    remove_rownames %>% tibble::column_to_rownames(var="isoform")
)


# phenotype
phenotype <- list(
  WholeTargeted = read.csv(paste0(dirnames$root_sfari, "/12_deseq2/WholeTargetedphenotype.csv"))
)

# subset phenotype
sphenotype <- list(
  TargetedGroup = phenotype$WholeTargeted %>% filter(grepl("Targeted",sample)) %>% mutate(time = age, col = paste0(sample,"_",group)),
  TargetedSex = phenotype$WholeTargeted %>% filter(grepl("Targeted",sample)) %>% mutate(time = age, col = paste0(sample,"_",sex)),
  WholeGroup = phenotype$WholeTargeted %>% filter(grepl("Whole",sample)) %>% mutate(time = age, col = paste0(sample,"_",group)),
  WholeSex = phenotype$WholeTargeted %>% filter(grepl("Whole",sample)) %>% mutate(time = age, col = paste0(sample,"_",sex))
)

# factors
factorsInput <- list(
  TargetedGroup = phenotype$WholeTargeted %>% filter((grepl("Targeted",sample))) %>% dplyr::select(sample, group) %>% magrittr::set_rownames(.$sample) %>%
    dplyr::rename("Replicate" = "group") %>%
    mutate(Replicate = ifelse(Replicate == "Prenatal",1,2)),
  TargetedSex = phenotype$WholeTargeted %>% filter((grepl("Targeted",sample))) %>% dplyr::select(sample, sex) %>% magrittr::set_rownames(.$sample) %>%
    dplyr::rename("Replicate" = "sex") %>%
    mutate(Replicate = ifelse(Replicate == "M",1,2)),
  WholeGroup = phenotype$WholeTargeted %>% filter((grepl("Whole",sample))) %>% dplyr::select(sample, group) %>% magrittr::set_rownames(.$sample) %>%
    dplyr::rename("Replicate" = "group") %>%
    mutate(Replicate = ifelse(Replicate == "Prenatal",1,2)),
  WholeSex = phenotype$WholeTargeted %>% filter((grepl("Whole",sample))) %>% dplyr::select(sample, sex) %>% magrittr::set_rownames(.$sample) %>%
    dplyr::rename("Replicate" = "sex") %>%
    mutate(Replicate = ifelse(Replicate == "M",1,2))
  
)
factorsInput <- lapply(factorsInput, function(x) x[order(x$Replicate),, drop = FALSE])

for(i in 1:4){
  write.table(factorsInput[[i]], paste0(dirnames$meta,"/",names(factorsInput)[[i]],"Factors.txt"), quote=F, sep="\t")
  write.csv(sphenotype[[i]], paste0(dirnames$meta,"/",names(sphenotype)[[i]],"Phenotype.csv"), quote=F)
}


## ---------- runDIU -----------------


resultsDIU <- list(
  ontTargGroup = runDIU(transMatrixRaw=Exp$ontTargTranRaw,transMatrix=Exp$ontTargTranNorm,classf=class.files$targ_SQ,
                        myfactors=factorsInput$targ_group,filteringType="FOLD",filterFC=2),
)

## ---------- write output -----------------

saveRDS(resultsDIU, file = paste0(dirnames$targ_output, "/resultsDIU.RDS"))

for(i in 1:3){
  write.csv(resultsDIU[[i]], paste0(dirnames$targ_output,"/",names(resultsDIU)[[i]],".csv"),row.names=F)
}