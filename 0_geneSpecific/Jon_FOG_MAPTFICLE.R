LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))


ficleDir <- "/lustre/projects/Research_Project-MRC148213/lsl693/SFARI/FICLE/SFARI"
classFilesName = "/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/SQANTI3_whole_RulesFilter_result_classification.targetgenes_counts_filtered.txt"
classFiles = SQANTI_class_preparation(classFilesName,"ns")
inputGtf = as.data.frame(rtracklayer::import(paste0(ficleDir, "/MAPT_sqanti3Filtered_classification.filtered_lite.gtf")))


A5A3 <- input_FICLE_splicing_results(ficleDir,"A5A3_tab")
pA5A3 <- plot_A5A3_Tgene(ficleDir, "MAPT") + theme(legend.position = "top")
pES <- plot_ES_Tgene(ficleDir,"MAPT")
pdenro <- plot_dendro_Tgene(ficleDir, "MAPT")

# track
exons <- inputGtf %>% dplyr::filter(type == "exon")
rescaled <- shorten_gaps(
  exons, 
  to_intron(exons, "transcript_id"), 
  group_var = "transcript_id"
)

track <- rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    fill = gene_id,
    y = transcript_id
  )) +
  geom_range(
  ) +
  geom_intron(
    data = rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200 
  ) + theme_classic() + scale_fill_manual(values = "blue") +
  theme(legend.position = "None") + labs(y = "Transcripts")


library("cowplot")
pdf("Mapt.pdf", width = 22, height = 15)
plot_grid(plot_grid(pdenro, pES,pA5A3, ncol = 1),track, nrow = 1) 
dev.off()
