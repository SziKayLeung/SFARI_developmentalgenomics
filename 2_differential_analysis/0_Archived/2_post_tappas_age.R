## ---------- Script -----------------
##
## Script name: 2_post_tappas_age
##
## Purpose of script: Generate plots post-tappAS of transcripts differentially expressed/spliced after accounting for age (fetal, adult as categories)
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## Performed Case-Control on tappAS (java on knight) with just age
##      Group --> Fetal (Control) and Adult (Case)
## Input files generated in 1_normalise_run_tappas.sh [A) Fetal vs Adult Brain] 
## Ouput files stored in sfari_differential.config.R


## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/"
source(paste0(SC_ROOT, "0_source_differential_functions.R"))
source(paste0(SC_ROOT, "sfari_differential.config.R"))


## ---------- Load tappAS files -----------------
loaded <- list(
  ont = input_tappasfiles(TAPPAS_INPUT_DIR$ont),
  ont_multi = input_tappasfiles(TAPPAS_INPUT_DIR$ont_multi)
)

# number of transcripts removed using tappAS filtering
filtered_p <- num_tappas_filter(loaded$ont$input_normalized_matrix, class.files)


## ---------- Annotate tappAS files -----------------
annotated <- list(
  ont = annotate_tappasfiles(class.files,loaded$ont$input_normalized_matrix,phenotype$ont),
  ont_multi = annotate_tappasfiles(class.files,loaded$ont_multi$input_normalized_matrix,phenotype$ont_multi)
)


## ---------- Differential Gene and Transcript results  -----------------

diff_gene  <- diff_results(loaded$ont,"gene",0.05)
diff_trans <- diff_results(loaded$ont,"transcript",0.05)

# WCPG 
WCPG_genes = c("CACNA1G","CUL1","GRIA3","GRIN2A","HERC1","RB1CC1", "SETD1A","SP4","TRIO","XPO7")
diff_trans[diff_trans$associated_gene %in% WCPG_genes,]

# Interesting genes from Rosie 
Rosie_genes = c("KIAA1109","SCN1A","SCN8A","KCNQ2","RBFOX3")
diff_trans[diff_trans$associated_gene %in% Rosie_genes,]


## ---------- Gene and Transcript expression plots  -----------------

# volcano plots
diff_gene_all <- plot_volcano(diff_gene)
diff_trans_all <- plot_volcano(diff_trans)


# gene expression plots
diff_gene_p <- list(
  top10 = generate_plots(unique(diff_gene_all$top10$associated_gene),annotated$ont,"Gene","ONT Transcript Expression"),
  wcpg = generate_plots(WCPG_genes,annotated$ont,"Gene","Iso-Seq Gene Expression")
)

# transcript expression plots
diff_trans_p <- list(
  top10 = generate_plots(unique(diff_trans_all$top10$associated_gene),annotated$ont,"Transcript","ONT Transcript Expression"),
  wcpg = generate_plots(WCPG_genes,annotated$ont,"Transcript","ONT Transcript Expression"),
  rosie = generate_plots(Rosie_genes,annotated$ont,"Transcript","ONT Transcript Expression")
)



## ---------- Differential Isoform Usage  -----------------

dev_genes <- c("VCAN","MAPT","RBFOX3")

IF_exp_rExp_p <- list()
for (i in 1:length(dev_genes)) {
  gene <- dev_genes[[i]]
  
  IF_exp_rExp_p[[i]] <- list(
    IF = IF_plot(gene,loaded$ont$gene_transcripts,loaded$ont$input_normalized_matrix, phenotype$ont, class.files),
    Exp = plot_trans_exp(gene,annotated$ont$Norm_transcounts,"All","ONT Expression"),
    rExp = plot_raw_expression_bygene(FL_reads, gene)
    )
}
names(IF_exp_rExp_p) <- dev_genes


pdf(paste0(output_dir,"/Gene_specific_plots.pdf"), onefile = TRUE, height = 12, width = 25)
plot_grid(plot_grid(IF_exp_rExp_p$VCAN$IF, IF_exp_rExp_p$VCAN$Exp), IF_exp_rExp_p$VCAN$rExp, nrow = 2)
plot_grid(plot_grid(IF_exp_rExp_p$MAPT$IF, IF_exp_rExp_p$MAPT$Exp), IF_exp_rExp_p$MAPT$rExp, nrow = 2)
plot_grid(plot_grid(IF_exp_rExp_p$RBFOX3$IF, IF_exp_rExp_p$RBFOX3$Exp), IF_exp_rExp_p$RBFOX3$rExp, nrow = 2)
dev.off()


## ---------- Output  -----------------

# Top 10
pdf(paste0(output_dir,"/Top10_diff_gene_plots.pdf"), onefile = TRUE)
for(i in diff_gene_p$top10){print(plot_grid(i))}
dev.off()

pdf(paste0(output_dir,"/Top10_diff_transcript_plots.pdf"), onefile = TRUE)
for(i in diff_trans_p$top10){print(plot_grid(i))}
dev.off()


# WCPG
pdf(paste0(output_dir,"/WCPG_diff_transcript_plots.pdf"), onefile = TRUE)
for(i in diff_trans_p$wcpg){print(plot_grid(i))}
dev.off()

pdf(paste0(output_dir,"/WCPG_raw_expression_plots.pdf"), onefile = TRUE, width = 30, heigh = 15)
for(i in WCPG_genes){print(plot_raw_expression_bygene(FL_reads, i))}
dev.off()


pdf(paste0(output_dir,"/Rosie_diff_transcript_plots.pdf"), onefile = TRUE, width = 30, heigh = 15)
for(i in diff_trans_p$rosie){print(plot_grid(i))}
dev.off()

pdf(paste0(output_dir,"/Rosie_raw_expression_plots.pdf"), onefile = TRUE, width = 30, heigh = 15)
for(i in Rosie_genes){print(plot_raw_expression_bygene(FL_reads, i))}
dev.off()
