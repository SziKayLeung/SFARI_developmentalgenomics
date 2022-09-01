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
## Performed Multiple Series Time-Course on tappAS (java on knight)
##      Time --> Fetal (0) and Adult (1)
##      Group --> Female (Control) and Male (Case)
## Input files generated in 1_normalise_run_tappas.sh  [B_Fetal vs Adult Brain (Sex differences)] 
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


## ---------- Transcript expression plots  -----------------
tgenes <- c("HHATL","CNDP1","OPALIN")
ttrans <- c("ENST00000441594.6","ENST00000358821.8","ENST00000371172.8")

'
Generate trajectory plots of gene (with all isoforms kept after normalisation): female vs male, fetal vs adult
Example: twocate_plot_transexp_overtime("HHATL",annotated$ont_multi$Norm_transcounts,"ONT_Expression")
twocate_plot_transexp_overtime <gene> <annotated_normalised_matrix> <plot_title>
'
diff_trans_overtime_p <- list(
  tgenes = generate_plots(tgenes,annotated$ont_multi,"2Cate_Transcript_overtime","ONT Expression")
)


'
Generate box-plots of specific transcript: female vs male, fetal vs adult
Example: twocate_plot_transexp("ENST00000441594.6",loaded$ont_multi$input_normalized_matrix,phenotype$ont_multi)
twocate_plot_transexp <transcript> <input_normalised_matrix> <phenotype_file>
'
diff_trans_p <- list(
  ttrans = generate_plots(ttrans,loaded$ont_multi$input_normalized_matrix,"2Cate_Transcript",phenotype$ont_multi) 
)


## ---------- Output  -----------------
# Top 3
pdf(paste0(output_dir,"/AgevsSex_diff_trans_plots.pdf"), onefile = TRUE, width = 17, height = 8)
for(i in (1:3)){
  print(plot_grid(diff_trans_overtime_p$tgenes[[i]], diff_trans_p$ttrans[[i]], rel_widths = c(0.7,0.3)))
}
dev.off()