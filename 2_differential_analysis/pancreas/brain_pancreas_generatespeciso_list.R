brain_pancreas_class.files <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/10_characterisation_SK/Brain_vs_Pancreas/2_sqanti3/merged_brain_pancreas_fetal_targeted_classification.txt", header = T, sep = "\t", as.is = T)

brain_pancreas_targetgenes <- brain_pancreas_class.files[brain_pancreas_class.files$associated_gene %in% targetgenes$V1,]
if(length(targetgenes$V1) != length(unique(brain_pancreas_targetgenes$associated_gene))){
  setdiff(targetgenes$V1, unique(brain_pancreas_targetgenes$associated_gene))
}

write.table(brain_pancreas_targetgenes$isoform, "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/9_tappAS_SK/1_Input/C_BrainvsPancreas/brain_pancreas_targetgenes_isoforms.txt", quote = F, row.names = F, col.names = F)

merged_brain_gff3 <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/10_characterisation_SK/Brain_vs_Pancreas/2_sqanti3/merged_brain_pancreas_fetal_targeted.gff3", as.is = T, sep = "\t")

merged_brain_gff3_target <- merged_brain_gff3 %>% filter(V1 %in% brain_pancreas_targetgenes$isoform)
write.table(merged_brain_gff3_target, "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted_targetgenes.gff3", sep = " ", quote = F, row.names = F, col.names = F)
