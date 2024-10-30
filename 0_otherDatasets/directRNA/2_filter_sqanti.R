#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: S.Leung (S.K.Leung@exeter.ac.uk)
## Previously R.Bamford processed direct RNA sequencing dataset using the Iso-Seq collapse pipeline for cDNA 
## given used as a validation dataset
## R script to do similar filtering of removing mono-exonic intergenic transcripts, mono-exonic transcripts associated to multi-exonic genes
## filter 2 reads, 2 samples
## --------------------------------


LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))

dirnames <- list(
  utils = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_utils/",
  directRNA = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/",
  demux = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/8_demux/"
)

class.names.files <- list(
  directRNA = paste0(dirnames$directRNA, "7_sqanti/merged_all_collapsed_corrected_filtered_RulesFilter_result_classification.txt")
)
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons != "Intergenic_1"))

# remove mono-exonic transcripts 
refExonNum <- read.csv(paste0(dirnames$utils,"gencode.v40.annotation.numExon.csv"))
exonicClass <- refExonNum %>% 
  group_by(Gene, GeneName, GeneType) %>% 
  dplyr::summarise(maxNumExon = max(NumExon)) %>% 
  mutate(monoExonic = ifelse(maxNumExon == 1, TRUE, FALSE))

refGeneType <- read.table(paste0(dirnames$utils,"gencode.v38.annotation.geneannotation.txt"), header = T)

# identify multi-exonic genes
multiExonicGenes <- unique(exonicClass[exonicClass$monoExonic == FALSE,][["GeneName"]])
# TRUE = monotranscript within multi-exonic gene
class.files <- lapply(class.files, function(x) x %>% mutate(monomulti = ifelse(exons == 1 & associated_gene %in% multiExonicGenes, TRUE,FALSE)))
# filter FALSE to retain multi-transcripts within multi-exonic gene, and mono-transcripts within mono-exonic gene
class.files <- lapply(class.files, function(x) x %>% filter(monomulti == FALSE))

# output classification file
write.table(class.files$directRNA, paste0(dirnames$directRNA, "9_sqanti_final/sqantifiltered_monoexonicfiltered_classification.txt"), sep = "\t", quote = F, row.names = F)

# filter 2 reads, 2 samples
demux_file <- fread(paste0(dirnames$demux, "merged_fl_count.csv"), data.table = FALSE)
rownames(demux_file) <- demux_file$id
demux_file <- demux_file %>% select(-id) %>% mutate(nsamples = rowSums(.!=0), nreads = rowSums(.)) 
class.files$directRNA <- merge(class.files$directRNA, demux_file, by.x = "isoform", by.y = "row.names", all.x = T)
class.files$directRNA <- class.files$directRNA %>% filter(nsamples >= 2, nreads >= 2)

# output classification file, further filtered by 2 reads 2 samples
write.table(class.files$directRNA, paste0(dirnames$directRNA, "9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt"), sep = "\t", row.names = F, quote = F)
write.table(class.files$directRNA$isoform, paste0(dirnames$directRNA, "9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_ID.txt"), quote = F, row.names = F, col.names = F)