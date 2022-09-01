## ---------- Script -----------------
##
## Script name: sfari_differential.config.R
##
## Purpose of script: Store the output files from tappAS for downsteram plots
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
## 
## ont = Fetal vs adult (age only), with case-control analysis --> plots in post_tappas_age.R
## ont_multi = Fetal vs adult, male vs female, with multiple series time-course --> plots in post_tappas_agevssex.R
##
## same classification file, given multiple analysis on the same data

ROOT_DIR = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/9_tappAS_SK/"
output_dir = paste0(ROOT_DIR, "01_figures_tables")

## ---------- TappAS input files -----------------

# output files from running tappAS
TAPPAS_INPUT_DIR = list(
  ont = paste0(ROOT_DIR,"2_Results/A_FetalvsAdult"),
  ont_multi = paste0(ROOT_DIR,"2_Results/B_FetalAdultvsSex")
)

# phenotype 
# column 1: sample; column 2: group
TAPPAS_PHENOTYPE = list(
  ont = paste0(ROOT_DIR, "1_Input/A_FetalvsAdult/phenotype.txt"),
  ont_multi = paste0(ROOT_DIR, "1_Input/B_FetalvsAdultvsSex/phenotype.txt")
)
phenotype <- lapply(TAPPAS_PHENOTYPE, function(x) read.table(x, header = T))


## ---------- SQANTI classification files -----------------

# Classification file
ISOSEQ_WKD_ROOT="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/8_SQANTI3_RB/"
class.names.files <- paste0(ISOSEQ_WKD_ROOT,"sqanti3Filtered2_classification.filtered_lite_classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"nstandard")


## ---------- Abundance -----------------

# raw FL reads
FL_reads <- read.table(paste0(ROOT_DIR, "1_Input/A_FetalvsAdult/SFARI_talon_expression.txt"), sep = "\t") 


