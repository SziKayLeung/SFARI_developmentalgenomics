#!/usr/bin/env Rscript
## ----------Script-----------------
## 
## SK.Leung: further filter by 2 reads, 2 samples (sqanti default expression file)
## merge sqanti filtered by mono-exonic transcripts with expression file (demux)
## ---------------------------------

suppressWarnings({
  suppressMessages(library("data.table"))
  suppressMessages(library("dplyr"))
})

# other function scripts
LOGEN <- c("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/")
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))

dir <- c("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_default/")

# SQANTI classification file
class.files.names <- c(paste0(dir, "WholeTargeted_collapsedAllChr_RulesFilter_result_classification.txt"))
class.files <- SQANTI_class_preparation(class.files.names, "ns")    

nrow(class.files)
# 7060521

dirnames <- list(
  output ="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/",
  readstat = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux/2_finalised_readstat",
  demux = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux/4_demux_merged"
)

# read in demux files across 24 chromosomes  
demux.names.files <- list.files(path = dirnames$demux, pattern = "fl_count.csv", full = TRUE)
demux.files <- lapply(demux.names.files, function(x) fread(x, data.table=FALSE))
merged.demux.files <- bind_rows(demux.files)

# filter demux files
filtered.merged.demux.files <- merged.demux.files %>% filter(nreads >= 2 & nsamples >= 2)
filtered.class.files <- class.files %>% filter(isoform %in% filtered.merged.demux.files$id)
nrow(filtered.class.files)
# 1881979

# final classification file 
final.class.files <- merge(filtered.class.files, filtered.merged.demux.files, by.x = "isoform", by.y = "id", all.x = TRUE)

# Loop through unique values in 'chrom'
for (i in unique(final.class.files$chrom)) {
  final.class.files <- final.class.files %>%
    mutate(associated_gene = if_else(
      str_detect(associated_gene, "novelGene") & chrom == i,
      str_replace(associated_gene, "novelGene_", paste0("novelGene_", i, "_")),
      associated_gene
    ))
}

write.table(final.class.files, paste0(dir,"sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt"), row.names = F, quote = F, 
            sep = "\t")