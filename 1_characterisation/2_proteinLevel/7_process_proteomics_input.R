#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: process output from long-read proteogenomics pipeline
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 

## ------------ packages ------------

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN, "merge_characterise_dataset/run_ggtranscript.R"))


## ------------ directory paths ------------

SC_ROOT <- "/lustre/projects/Research_Project-MRC148213/lsl693/SFARI_developmentalgenomics/1_characterisation"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
dirnames <- list(
  sqanti = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/",
  proteomics = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics/",
  utils = paste0(root_sfari,"/0_utils/")
)


## ------------ input ------------

input <- list(
  # transcript SQANTI3 classification file
  t.class.files = SQANTI_class_preparation(paste0(dirnames$sqanti,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic_modFL.txt"),"ns"),
  # protein unfiltered long-read proteogenomics classification file
  #pUnfiltered.class.files = pSQANTI_class_preparation(paste0(dirnames$proteomics,"7_classified_protein/ADBDR_unfiltered.protein_classification.tsv")),
  # protein filtered long-read proteogeonomics classification file
  #Filtered.class.files = pSQANTI_class_preparation(paste0(dirnames$root,"7_classified_protein/ADBDR.sqanti_protein_classification.tsv")),
  # list of transcripts collapsed by protein reading frame
  t2p.collapse = fread(paste0(dirnames$proteomics,"6_refined_database/Whole_orf_refined.tsv")),
  # gtf of aligned peptide ORF
  #peptide_orf = paste0(dirnames$root, "5_calledOrfs/ADBDR.gtf"),
  # cpat
  cpat_output = fread(paste0(dirnames$proteomics,"5_calledOrfs/Whole.ORF_prob.best.tsv")),
  no_cpat_orf = fread(paste0(dirnames$proteomics, "5_calledOrfs/Whole.no_ORF.txt"),header = F)
)


# number of RNA transcripts collapsed by protein sequence
input$t2p.collapse <- input$t2p.collapse %>% mutate(numtxCollapsed = count.fields(textConnection(pb_accs), sep = "|"))
char <- strsplit(as.character(input$t2p.collapse$pb_accs), '|', fixed = T)
t2p.collapse.dissected <- data.frame(pb_accs=unlist(char), base_acc=rep(input$t2p.collapse$base_acc, sapply(char, FUN=length)))
input$t2p.collapse <- merge(t2p.collapse.dissected,input$t2p.collapse[,c("base_acc","numtxCollapsed")])
input$t.class.files <- merge(input$t.class.files, input$t2p.collapse, by.x = "isoform", by.y = "pb_accs")


## ------------ output ------------ 

proteinInput <- input
save(proteinInput, file = paste0(dirnames$utils,"/proteinInputWhole.RData"))

# aggregate sum by same peptide sequence
pFL <- input$t.class.files %>% dplyr::select(base_acc, contains("FL."))
pFLsum <- aggregate(. ~ base_acc, pFL, sum)