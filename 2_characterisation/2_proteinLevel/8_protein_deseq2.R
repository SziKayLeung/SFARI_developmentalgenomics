## ---------- Script -----------------
##
## Purpose: perform differential analysis on AD-BDR Iso-Seq targeted datasets using linear regression
## Transcript level separate analysis
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
##
## ---------- Notes ------------------
## two ONT batches: Batch 1 - Nov 2022, Batch 2 - March 2023 


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))


## ---------- source functions -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))

label_group <- function(genotype){
  if(genotype %in% c("Case","CASE")){group = "AD"}else{
    if(genotype %in% c("Control","CONTROL")){group = "Control"}}
  return(group)
}

## ---------- input -----------------

# directory names
SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation"
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  sqanti = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/",
  proteomics = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/6a_longReadProteogenomics/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics"
)


# classfiles
class.files = SQANTI_class_preparation(paste0(dirnames$sqanti,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic_modFL.txt"),"ns")

# proteomics input 
load(file = paste0(SC_ROOT,"/proteinInputWhole.RData"))

input <- list()
# phenotype
input$wholePhenotype <- fread(paste0(root_dir, "RBFetal/00_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% 
  mutate(group = factor(group, levels = c("Prenatal","Postnatal")), col = paste0(sample,"_",group)) %>% 
  filter(grepl("Whole",sample))
input$classfiles <- class.files %>% select(isoform, associated_gene, associated_transcript, structural_category, subcategory, contains("Whole"))
row.names(input$classfiles) <- input$classfiles$isoform


# phenotype
phenoFiles <- list(
  wholeGroup = input$wholePhenotype,
  wholeSex = input$wholePhenotype %>% mutate(group = sex)
)

### protein 
# aggregate sum by same peptide sequence
counts <- class.files %>% select(isoform, contains("Whole"))
pFL <- proteinInput$t.class.files %>% dplyr::select(base_acc, contains("Whole")) 
pFLsum <- aggregate(. ~ base_acc, pFL, sum)

# datawrangle for input to run_DESeq2()
expressionFiles <- list(
  iso = input$classfiles %>% dplyr::select(starts_with("Whole")),
  isoProtein = pFLsum %>% tibble::column_to_rownames(., var = "base_acc")
)
row.names(expressionFiles$iso) <- input$classfiles$isoform


## ---------- ONT: Creating DESeq2 object and analysis -----------------

ResTran <- list(
  #tWald = run_DESeq2(test="Wald",expressionFiles$iso,phenoFiles$isoGrp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  pWald = run_DESeq2(test="Wald",expressionFiles$isoProtein,phenoFiles$wholeGroup,threshold=10,controlname="Prenatal",design="case_control",groupvar="factor")
)

annoResTran <- list(
  #Wald = anno_DESeq2(ResTran$tWald,input$classfiles,phenoFiles$isoGrp,controlname="Control",level="transcript",sig=0.05),
  pWald = anno_DESeq2(ResTran$pWald,input$classfiles,phenoFiles$wholeGroup,controlname="Prenatal",level="transcript",sig=0.05)
)


## ---------- Output -----------------

saveRDS(annoResTran, file = paste0(dirnames$output, "/IsoProtein_Whole_DESeq2TranscriptLevel.RDS"))
