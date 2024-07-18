library("data.table")
library("ggplot2")
library("ggtranscript")
library("cowplot")

dat <- fread("/lustre/projects/Research_Project-MRC148213/lsl693/SFARI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt")
counts <- fread("/lustre/projects/Research_Project-MRC148213/lsl693/SFARI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic_counts.txt")

gtf <- as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC148213/lsl693/SFARI/SNCA_WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf"))

refgtf <- as.data.frame(rtracklayer::import(("/lustre/projects/Research_Project-MRC148213/lsl693/reference/annotation/SNCA_gencodeV40.gtf"))) %>% filter(gene_name == "SNCA")


gtfmerged <- rbind(gtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                    refgtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")]) 

SNCA <- dat[dat$associated_gene == "SNCA",]
SNCA <- counts[counts$associated_gene == "SNCA",]

nrow(SNCA)
# 60 SNCA transcripts differentially expressed associated with developmental age

phenotype <- read.csv("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/WholeTargetedphenotype.csv")

snca_DTE <- c("ONT4_6374_2774","ONT4_6374_2772","ONT4_6374_1807","ONT4_6374_2773","ONT4_6374_1811",
              "ONT4_6374_2761", "ONT4_6374_2815","ONT4_6374_1853")

pLifeTime <- list()
for(x in snca_DTE){
  print(x)
  pLifeTime[[x]] = plot_trans_exp_lifetime(x,class.files$glob_SQ,Exp$whole_group)
}
plot_grid(plotlist = pLifeTime)

plot_DTE <- function(transcripts){
  print(transcripts)
  dat <- SNCA %>% filter(isoform == transcripts) %>% select(contains("Whole")) %>% reshape2::melt() %>% 
  merge(., phenotype, by.x = "variable", by.y = "sample") %>% 
  mutate(group = factor(group, levels = c("Prenatal","Postnatal"))) 
  
  p <- ggplot(dat, aes(x = group, y = value)) + geom_boxplot() +
  labs(x = "", y = "FL reads", title = transcripts) + theme_classic()
  
  return(p)
}


scna_exon <- gtfmerged %>% filter(type == "exon") %>% 
  filter(transcript_id %in% c(snca_DTE,"ENST00000394991.8"))

rescaled <- shorten_gaps(
  scna_exon, 
  to_intron(scna_exon, "transcript_id"), 
  group_var = "transcript_id"
) %>% merge(., SNCA[,c("isoform","structural_category")],by.x= "transcript_id", by.y = "isoform", all.x = T)


B <- rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(aes(fill = structural_category)
  ) +
  geom_intron(
    data = rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200
  ) + theme_classic()


pDTEs <- lapply(snca_DTE, function(x) plot_DTE(x))
A <- plot_grid(plotlist = pDTEs)

plot_grid(A, B)
