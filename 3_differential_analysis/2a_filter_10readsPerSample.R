#!/usr/bin/env Rscript

# Szi Kay Leung
# 12.07.2024: subset expression file from whole and targeted dataset to ensure only isoforms with 10 reads in every sample is kept

suppressMessages(library("dplyr"))

# input
message("Read demux file")
demux <- data.table::fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/10_deseq/WholeTargeted_demux.csv", data.table = F)

message("Read classification file")
inputClassfile <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples_monomultirem.txt", sep = "\t", header = T)

message("Keep only final filtered isoforms")
filtered_demux <- demux[demux$isoform %in% inputClassfile$isoform, ]

# regenerate as matrix file, rownames = isoform
row.names(filtered_demux) <- filtered_demux$isoform
filtered_demux <- filtered_demux %>% select(-isoform)

# keep whole samples only 
whole_filtered_demux <- filtered_demux %>% dplyr::select(contains("Whole"))
whole_filtered_demux <- as.matrix(whole_filtered_demux)

# function to identify rows with minimum 10 reads across each column
message("Keeping isoforms with minimum 10 reads across each sample")
threshold = 10
row_check <- function(row, threshold) {
  all(row >= threshold)
}

# Apply the function to each row and keep rows where the condition is TRUE
expressionfiltered_demux <- whole_filtered_demux[apply(whole_filtered_demux, 1, row_check, threshold), ]

# output
write.csv(expressionfiltered_demux, "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/10_deseq/Whole_demux_10readsperSample.csv", row.names = T, quote = F)
