#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Rebuttal figures
##
## ---------------------------------

suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
phenotype <- data.table::fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)


## ------- rarefaction curves ------- 

# input files
#input_rarefaction_dir="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/0_rarefactionSKL" 
input_rarefaction_dir=paste0(root_dir,"rarefaction")
input_rarefaction_genes <- list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbgene\\.min_fl_2\\.txt", full.names = T)
input_rarefaction_transcripts <- list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbid\\.min_fl_2\\.txt", full.names = T)

rareFactionSamples <- intersect(
  word(list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbgene\\.min_fl_2\\.txt"),c(1),sep=fixed("_")),
  word(list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbid\\.min_fl_2\\.txt"),c(1),sep=fixed("_"))
)
manifest[manifest$ID %in% rareFactionSamples,]

rarefaction_genes <- lapply(input_rarefaction_genes, function(x) data.table::fread(x))
rarefaction_transcripts <- lapply(input_rarefaction_transcripts, function(x) data.table::fread(x))
names(rarefaction_genes) <- word(list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbgene\\.min_fl_2\\.txt"),c(1),sep = fixed("_"))
names(rarefaction_transcripts) <- word(list.files(path = input_rarefaction_dir, pattern = "Whole.*by_pbid\\.min_fl_2\\.txt"),c(1),sep = fixed("_"))
merged_genes <- bind_rows(rarefaction_genes, .id = "sample") 
merged_transcripts <- bind_rows(rarefaction_transcripts, .id = "sample") 

# plot distribution of genes
merged_genes <- merge(merged_genes, phenotype, by = "sample", all.x = T) %>% filter(sample %in% rareFactionSamples)
merged_transcripts <- merge(merged_transcripts, phenotype, by = "sample", all.x = T) %>% filter(sample %in% rareFactionSamples)

pRarefaction1 <- ggplot(merged_genes , aes(x = size, y = mean, colour = sample, linetype = group)) + 
  geom_line(size = 1.5) + 
  labs(x ="Number of subsampled reads (K)", y = "Number of genes", linetype = NULL) + 
  theme_classic() + theme(legend.position = "bottom") +
  guides(colour = "none")

pRarefaction2 <- ggplot(merged_transcripts , aes(x = size, y = mean, colour = sample, linetype = group)) + 
  geom_line(size = 1.5) + 
  labs(x ="Number of subsampled reads (K)", y = "Number of transcripts", linetype = NULL) + 
  theme_classic() + theme(legend.position = "bottom") +
  guides(colour = "none")

pdf(paste0(output_dir,"rarefaction.pdf"),width = 15, height = 8)
plot_grid(pRarefaction1, pRarefaction2, nrow = 1, labels = c("A","B"))
dev.off()

## ------ raw reads ---------

manifest <- read.csv("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/WholeTargetedphenotype_manifest.csv")
clusterReport <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/4_transcriptClean/cluster_report_counts.txt")
clusterReport <- clusterReport %>% mutate(sample = word(V2,c(1),sep=fixed("."))) 
merge(clusterReport, manifest, by.x = "sample", by.y = "ID")



