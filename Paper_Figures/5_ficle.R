LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))


ficleDir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp"

pGRIK2 <- list(
  overview = plot_dendro_Tgene(ficleDir, "GRIK2"),
  A5A3 = plot_A5A3_Tgene(ficleDir, "GRIK2"),
  ES = plot_ES_Tgene(ficleDir,"GRIK2",class.files$glob_targ_SQ)[[1]],
  tracks = p
)

pTRIO <- list(
  overview = plot_dendro_Tgene(ficleDir, "TRIO"),
  A5A3 = plot_A5A3_Tgene(ficleDir, "TRIO"),
  ES = plot_ES_Tgene(ficleDir,"TRIO",class.files$glob_targ_SQ)[[1]],
  tracks = p2
)

SFARI_genes <- c("SHANK3","MECP2","CHD8","CNATNAP2","PTEN","SCN1A","SYNGAP1","ADNP","ARD1B","DYRK1A")
class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene %in% SFARI_genes,"isoform"]

WholeDESeqSig$sex %>% filter(!chrom %in% c("chrX","chrY")) %>% 
  left_join(., class.files$glob_targ_SQ[,c("isoform","associated_gene")], by = "isoform") %>% 
  mutate(SFARI_1_2 = ifelse(associated_gene %in% SFARI_CLASS_1_2,TRUE,FALSE),
         SFARI_1_2_S = ifelse(associated_gene %in% SFARI_CLASS_1_2_S ,TRUE,FALSE))

plot_grid(plot_trans_exp_individual("ONT5_6069_349",class.files$glob_targ_SQ,Exp$whole_sex,"sex") + facet_grid(~group),
          plot_trans_exp_lifetime("ONT5_6069_349",class.files$glob_targ_SQ,Exp$whole_sex)
)

plot_grid(
  plot_grid(pGRIK2[[1]],pGRIK2[[3]],pGRIK2[[2]] + theme(legend.position="top"),ncol=1),
  pGRIK2[[4]]
)

plot_grid(
  plot_grid(pTRIO[[1]],pTRIO[[3]] + scale_x_discrete(breaks = seq(10, 60, by = 10)),pTRIO[[2]] + 
              scale_x_discrete(breaks = seq(10, 60, by = 10)) + 
              theme(legend.position="top"),ncol=1),
  pTRIO[[4]], rel_widths = c(0.5,0.5)
)
