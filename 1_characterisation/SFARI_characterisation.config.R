## ---------- Script -----------------
##
## Script name: SFARI_characterisation.config.R
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
## 
## 
## 
## 
## 

ROOT_DIR = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/"
PANCREAS_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0055_20220623_10661/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined"
output_dir = paste0(ROOT_DIR, "9_tappAS_SK/01_figures_tables")


## ---------- Packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/characterise/sqanti_general.R")


## ---------- SQANTI classification files -----------------

class.files.names <- list(
  brain = paste0(ROOT_DIR,"8_SQANTI3_RB/sqanti3Filtered2_classification.filtered_lite_classification.txt"),
  pancreas = paste0(PANCREAS_DIR,"Targeted_AM_SQANTI3/Targeted_AM_SQANTI3_classification.filtered_lite_classification.txt"),
  brain_pancreas = paste0(ROOT_DIR,"10_characterisation_SK/Brain_vs_Pancreas/2_sqanti3/SFARI.annotated_classification.txt")
) 

class.files <- lapply(class.files.names, function(x) SQANTI_class_preparation(x,"nstandard"))


## ---------- CPAT related plots -----------------

CPAT.files.names <- list(
  brain = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/scripts/CPAT/10380_whole.ORF_prob.best.tsv"
)

CPAT <- lapply(CPAT.files.names, function(x) read.table(x, sep = "\t", header = T))



## ---------- Gffcompare: Brain vs Pancreas -----------------

brainvspancreas <- list(
  # generated from brain_pancreas_mergedexpression.R
  gffcomp_match_output = read.table(paste0(ROOT_DIR, "9_tappAS_SK/1_Input/C_BrainvsPancreas/brain_pancreas_common_genes_isoforms.txt"), header = T),
  gffcomp_expression = read.table(paste0(ROOT_DIR, "9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_expression.txt"), header = T) %>% 
    tibble::rownames_to_column(., "isoform") ,
  phenotype = read.table(paste0(ROOT_DIR, "9_tappAS_SK/1_Input/C_BrainvsPancreas/brain_pancreas_phenotype.txt"), header = T)
)

# reannotate brain and pancreas class files to include results from gffcompare 
class.files$brain <- merge(class.files$brain, brainvspancreas$gffcomp_match_output[,c("brain_isoform","brain_pancreas_common_isoform")], 
                           by.x = "isoform", by.y = "brain_isoform", all.x = T) %>% 
  mutate(common_unique_isoform = ifelse(is.na(brain_pancreas_common_isoform),isoform, as.character(brain_pancreas_common_isoform)))

class.files$pancreas <- merge(class.files$pancreas, brainvspancreas$gffcomp_match_output[,c("pancreas_isoform","brain_pancreas_common_isoform")], 
                           by.x = "isoform", by.y = "pancreas_isoform", all.x = T) %>%
  mutate(common_unique_isoform = ifelse(is.na(brain_pancreas_common_isoform),isoform, as.character(brain_pancreas_common_isoform)))

