#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Aim: Extract TMEM106B isoforms and respective ONT raw reads from WholeTargeted dataset
## Colloboration: Alex Salazar, Henne Holstege
## --------------------------------

library("data.table")
library("dplyr")

# load the classification file
root_rb_dir <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/"
Tmem_dir = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_colloborations/TMEM106B/"
dirnames <- list(wholetarg_SQ = paste0(root_rb_dir,"6_sqanti3/"))
load(file = paste0(dirnames$wholetarg_SQ,"all_filtered_classification_2reads2samples_noMonoIntergenicAll.RData"))
class.files$glob_targ_SQ_counts = fread(paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic_counts.txt"), data.table = F)

# extract the TMEM106B isoforms
TMEM106BIsoform <- class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene == "TMEM106B","isoform"]


# extract TMEM106B classification related information
TMEM106B.class.files <- class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene == "TMEM106B",
                                                 c("isoform","associated_gene","associated_transcript","structural_category")] 

TMEM106B.counts <- class.files$glob_targ_SQ_counts %>% filter(isoform %in% TMEM106BIsoform) %>% select(isoform, nreads,nsamples)
TMEM106B.class.files <- merge(TMEM106B.class.files, TMEM106B.counts, by = "isoform")

# output
write.table(TMEM106BIsoform, paste0(Tmem_dir, "TMEM106B_isoforms.txt"), col.names = F, row.names = F, quote = F)
write.table(TMEM106B.class.files,paste0(Tmem_dir, "TMEM106B_classification.txt"), col.names = F, row.names = F, quote = F)

## bash 
#sfari=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI
#grep -w -f ${sfari}/6_sqanti/sqanti/TMEM106B_isoforms.txt ${sfari}/5_isoseq/WholeTargeted/WholeTargeted_cleaned_aligned_merged_collapsed_chr7.read_stat.renamed.txt > ~/TMEM106B_WholeTargeted.read_stat.txt

