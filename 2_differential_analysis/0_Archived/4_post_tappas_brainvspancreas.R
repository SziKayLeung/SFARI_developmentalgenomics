## ---------- Script -----------------
##
## Script name: 3_post_tapas_agevssex.R
##
## Purpose of script: Generate plots post-tappAS of transcripts differentially expressed/spliced after accounting for age and sex
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
## 
## Performed Case control analysis on tappAS (java on knight)
##      Group --> Brain (Control) and Pancreas (Case)
## Targeted fetal datasets      
## Input files generated in 1_normalise_run_tappas.sh  [B_Fetal vs Adult Brain (Sex differences)] 
## Ouput files stored in sfari_differential.config.R


## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/"
source(paste0(SC_ROOT, "0_source_differential_functions.R"))
source(paste0(SC_ROOT, "sfari_differential.config.R"))


## ---------- Load tappAS files -----------------
loaded <- list(
  ont = input_tappasfiles(TAPPAS_INPUT_DIR$ont),
  ont_multi = input_tappasfiles(TAPPAS_INPUT_DIR$ont_multi),
  ont_tissue = input_tappasfiles(TAPPAS_INPUT_DIR$ont_tissue)
)


## ---------- Annotate tappAS files -----------------
annotated <- list(
  ont = annotate_tappasfiles(class.files$brain_whole,loaded$ont$input_normalized_matrix,phenotype$ont),
  ont_multi = annotate_tappasfiles(class.files$brain_whole,loaded$ont_multi$input_normalized_matrix,phenotype$ont_multi),
  ont_tissue = annotate_tappasfiles(class.files$brain_pancreas_targeted,loaded$ont_tissue$input_normalized_matrix,phenotype$ont_tissue)
)


## ---------- Transcript expression plots  -----------------
# number of upregulated vs downregulated transcript expression 
loaded$ont_tissue$results_trans <- loaded$ont_tissue$results_trans %>% mutate(direction = ifelse(log2FC < 0, "Downregulated","Upregulated"),
                                                                              revprob = 1-prob)
loaded$ont_tissue$results_trans %>% group_by(direction) %>% tally()

# top10 differentially expressed transcripts
# MAPT 
MAPT_isoforms <- class.files$brain_pancreas_targeted %>% filter(associated_gene == "MAPT") %>% .[,"isoform"]
MAPT_isoforms_diff <- loaded$ont_tissue$results_trans %>% filter(transcript %in% MAPT_isoforms) %>% arrange(log2FC)

diff_trans <- list(
  top10 = loaded$ont_tissue$results_trans %>% arrange(revprob) %>% .[1:10,"transcript"],
  top10_up = loaded$ont_tissue$results_trans %>% filter(direction == "Upregulated") %>% arrange(-log2FC) %>% .[1:10,"transcript"],
  top10_down = loaded$ont_tissue$results_trans %>% filter(direction == "Downregulated") %>% arrange(log2FC) %>% .[1:10,"transcript"],
  pancreas_only = loaded$ont_tissue$results_trans %>% filter(X1_mean < 5 & X2_mean > 10) %>% arrange(log2FC) %>% .[1:10,"transcript"],
  brain_only = loaded$ont_tissue$results_trans %>% filter(X2_mean < 5 & X1_mean > 10) %>% arrange(-log2FC) %>% .[1:10,"transcript"],
  brain_pancreas = loaded$ont_tissue$results_trans %>% filter(X1_mean > 50 & X2_mean > 50) %>% arrange(revprob) %>% .[1:10,"transcript"],
  mapt =  loaded$ont_tissue$results_trans %>% filter(transcript %in% MAPT_isoforms) %>% arrange(log2FC) %>% .[,"transcript"]
)

# transcript expression plots
diff_trans_p <- lapply(diff_trans, function(x) generate_plots(x, annotated$ont_tissue, "Per_Transcript","ONT Transcript Expression"))

# output plots
p1 = plot_grid(plotlist = diff_trans_p$top10)
p2 = plot_grid(plotlist = diff_trans_p$top10_up)
p3 = plot_grid(plotlist = diff_trans_p$top10_down)
p4 = plot_grid(plotlist = diff_trans_p$pancreas_only)
p5 = plot_grid(plotlist = diff_trans_p$brain_only)
p6 = plot_grid(plotlist = diff_trans_p$brain_pancreas)
p7 = plot_grid(plotlist = diff_trans_p$mapt)


## ---------- Transcript usage plots  -----------------

# genes with differential isoform usage with major isoforomo switch, arranged by qValue and total change
diff_iu <- list(
  PodChange = loaded$ont_tissue$results_DIU %>% filter(podiumChange == "TRUE") %>% arrange(qValue, -totalChange) %>% .[1:10,"gene"]
)

# transcript usage plots
diff_iu_p <- lapply(diff_iu$PodChange, function(x) generate_plots(x,loaded$ont_tissue,"usage",phenotype$ont_tissue, class.files$brain_pancreas_targeted))

## Examples
#IF_plot("GPM6A",loaded$ont_tissue$gene_transcripts,loaded$ont_tissue$input_normalized_matrix, phenotype$ont_tissue, class.files$brain_pancreas_targeted)
#IF_plot("MAGI2",loaded$ont_tissue$gene_transcripts,loaded$ont_tissue$input_normalized_matrix, phenotype$ont_tissue, class.files$brain_pancreas_targeted)


## ---------- Examples for BaseScope  -----------------

# examples of highly expressed genes with DIU
#highly_expressed_genes = brainvspancreas$mean_gene_exp %>% filter(log2Brain > 10 & log2Pancreas > 10)
#highly_expressed_DIU_genes = intersect(highly_expressed_genes$associated_gene, loaded$ont_tissue$results_DIU$gene)
#loaded$ont_tissue$results_DIU %>% filter(gene %in% highly_expressed_DIU_genes) %>% arrange(-totalChange)

#diff_IU <- lapply(highly_expressed_DIU_genes, function(x)
#                 IF_plot(x,loaded$ont_tissue$gene_transcripts,loaded$ont_tissue$input_normalized_matrix, phenotype$ont_tissue, class.files$pancreas)
#)

#plot_grid(plotlist = diff_IU)


pdf(paste0(output_dir,"/SFARI_brainvspancreas_differential.pdf"), width = 14, height = 10)
p1
p2
p3
p6
dev.off()

#pdf(paste0(output_dir,"/SFARI_brainvspancreas_differential_usage.pdf"), onefile = TRUE)
pdf(SFARI_brainvspancreas_differential_usage.pdf, onefile = TRUE)
for(i in 1:length(diff_iu_p)){
  print(plot_grid(diff_iu_p[[i]][[1]]))
}
dev.off()


# MAPT 
MAPT_isoforms <- class.files$brain_pancreas_targeted %>% filter(associated_gene == "MAPT") %>% .[,"isoform"]
MAPT_isoforms_diff <- loaded$ont_tissue$results_trans %>% filter(transcript %in% MAPT_isoforms) %>% arrange(log2FC)
loaded$ont_tissue$results_gene %>% filter(gene == "MAPT")

pdf(paste0(output_dir,"/SFARI_brainvspancreas_mapt_differential.pdf"), width = 14, height = 10)
plot_gene_exp("MAPT",annotated$ont_tissue$GeneExp,annotated$ont_tissue$Norm_transcounts,"Gene Expression")
p7
IF_plot("MAPT",loaded$ont_tissue$gene_transcripts, loaded$ont_tissue$input_normalized_matrix,phenotype$ont_tissue, class.files$brain_pancreas_targeted)
dev.off()
