## ---------- Script -----------------
##
## Script name: 
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

## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/"
source(paste0(SC_ROOT, "0_source_differential_functions.R"))
source(paste0(SC_ROOT, "sfari_differential.config.R"))

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/comp_characterise.R")
anno_dir = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/10_characterisation_SK"
MAPT_iso_class = read.csv(paste0(anno_dir, "/MAPT/Stats/MAPT_further_classifications.csv"))

## ---------- Load tappAS files -----------------
loaded <- list(
  ont = input_tappasfiles(TAPPAS_INPUT_DIR$ont)
)

# number of transcripts removed using tappAS filtering
filtered_p <- num_tappas_filter(loaded$ont$input_normalized_matrix, class.files)


## ---------- Annotate tappAS files -----------------
annotated <- list(
  ont = annotate_tappasfiles(class.files,loaded$ont$input_normalized_matrix,phenotype$ont)
)

pHeat <- list(
  mapt = draw_heatmap_gene("MAPT", class.files, annotated$ont$Norm_transcounts)
)
names(pHeat) = "MAPT"

plot_grid(pHeat$MAPT$gtable)
dendro_plot("MAPT")
pMAPT_class <- list(plot_mapt_classification("Human","sfari_ont"),plot_mapt_classification_bygroup())

plot_grid(pMAPT_class[[2]],pMAPT_class[[1]][[2]],rel_widths = c(0.6,0.4))
