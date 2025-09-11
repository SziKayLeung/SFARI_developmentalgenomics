#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Tabulate the number of protein isoforms from the list of selected target genes (n = 241) for Alessia Fgx Grant
## --------------------------------


## ---------- source config -----------------

SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))


## ---------- filter and tabulate -----------------

df <- class.files$protein_filtered %>% filter(associated_gene %in% selectedTargetGenes) %>%
  group_by(associated_gene, pr_splice_cat) %>% tally() %>%
  filter(pr_splice_cat %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>%
  spread(., key = pr_splice_cat, value = n) %>% 
  as.data.frame() %>%
  replace(is.na(.), 0)
colnames(df) <- c("associated_gene","pFSM","pISM","pNIC","pNNC")

write.csv(df,"C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/PaperZenodo/proteomics/targetGenes_numProteinsClassification.csv", row.names = F)
