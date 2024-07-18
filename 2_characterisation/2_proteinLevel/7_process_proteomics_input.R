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
  t.class.files = SQANTI_class_preparation(paste0(dirnames$sqanti,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic_modFL.txt"),"ns"),
  cpat = fread(paste0(dirnames$proteomics,"5_calledOrfs/Whole.ORF_prob.best.tsv"), data.table = F, header = T),
  cpat_best = fread(paste0(dirnames$proteomics,"5_calledOrfs/Whole_best_orf.tsv"), data.table = F, header = T),
  mapped = fread(paste0(dirnames$proteomics,"5_calledOrfs/all_orfs_mapped.tsv"), data.table = F, header = T),
  no_orf = fread(paste0(dirnames$proteomics, "5_calledOrfs/Whole.no_ORF.txt"), header = F),
  t2p.collapse = fread(paste0(dirnames$proteomics,"6_refined_database/Whole_orf_refined.tsv"), header = T)
)


# number of RNA transcripts collapsed by protein sequence
input$t2p.collapse.refined <- input$t2p.collapse %>% mutate(numtxCollapsed = count.fields(textConnection(pb_accs), sep = "|"))
char <- strsplit(as.character(input$t2p.collapse.refined$pb_accs), '|', fixed = T)
t2p.collapse.dissected <- data.frame(pb_accs=unlist(char), base_acc=rep(input$t2p.collapse.refined$base_acc, sapply(char, FUN=length)))
input$t2p.collapse.refined <- merge(t2p.collapse.dissected,input$t2p.collapse.refined[,c("base_acc","numtxCollapsed")])

## re-determine representative colalsped isoform: using ONT abundance (sum across all samples) rather than arbitrary (G.Shenkyman pipeline)
# take the ONT_sum read counts from the classification file
# max = grouping by the base_acc (i.e. the previously selected isoform), select the rows with the maximum ONT FL reads
# create an index to remap and create a "corrected_acc" column with the corresponding isoform that has the highest number of ONT FL reads
input$t2p.collapse.refined <- merge(input$t2p.collapse.refined,input$t.class.files[,c("isoform","FL")],by.x = "pb_accs", by.y = "isoform", all.x = TRUE)
max = input$t2p.collapse.refined %>% group_by(base_acc) %>% filter(FL == max(FL))
idx <- match(input$t2p.collapse.refined$base_acc, max$base_acc)
input$t2p.collapse.refined = transform(input$t2p.collapse.refined , corrected_acc = ifelse(!is.na(idx), as.character(max$pb_accs[idx]), base_acc))

## include in the original classification file the collapsed PB.ID 
input$t.class.files <- merge(input$t.class.files, input$t2p.collapse.refined[,c("pb_accs","numtxCollapsed","base_acc","corrected_acc")], by.x = "isoform", by.y = "pb_accs", all.x = TRUE)
 

## ------------ output ------------ 

proteinInput <- input
save(proteinInput, file = paste0(dirnames$utils,"/proteinInputWhole.RData"))

# aggregate sum by same peptide sequence
pFL <- input$t.class.files %>% dplyr::select(base_acc, contains("FL."))
pFLsum <- aggregate(. ~ base_acc, pFL, sum)