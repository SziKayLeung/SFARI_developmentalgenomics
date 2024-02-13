
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs/SuppTables")
dat <- data.frame()

for(i in 1:length(TargetGene)){
  gene = as.character(TargetGene[i])
  dat[i,1] <- gene
  dat[i,2] <- ifelse(nrow(WholeDESeqGeneSig$age[WholeDESeqGeneSig$age$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,3] <- ifelse(nrow(TargetedDESeqGeneSig$age[TargetedDESeqGeneSig$age$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,4] <- ifelse(nrow(WholeDESeqGeneSig$sex[WholeDESeqGeneSig$sex$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,5] <- ifelse(nrow(TargetedDESeqGeneSig$sex[TargetedDESeqGeneSig$sex$associated_gene == gene,]) == 1,TRUE,FALSE)
  dat[i,6] <- nrow(WholeDESeqSig$age %>% filter(associated_gene == gene))
  dat[i,7] <- nrow(TargetedDESeqSig$age %>% filter(associated_gene == gene))
  dat[i,8] <- nrow(WholeDESeqSig$sex %>% filter(associated_gene == gene))
  dat[i,9] <- nrow(TargetedDESeqSig$sex %>% filter(associated_gene == gene))
}
colnames(dat) <- c("TargetGene","WholeDGEGroup","TargetedDGEGroup","WholeDGESex","TargetedDGESex","NumWholeDTEGroup","NumTargetedDTEGroup","NumWholeDTESex","NumTargetedDTESex")

dat <- dat %>% mutate(DGEGroup = ifelse(WholeDGEGroup == TargetedDGEGroup, TRUE,FALSE),
               DGESex = ifelse(WholeDGESex == TargetedDGESex, TRUE,FALSE))


write.csv(dat,paste0(output_dir,"/TargetGenesDEAnalysis.csv"))
plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$whole_group,gene="MECP2")
plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$targeted_group,gene="MECP2")

plot_grid(
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$whole_group,"group","MECP2") + ylim(2,3),
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$targeted_group,"group","MECP2") + ylim(2,3)
)

plot_grid(
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$whole_group,"group","EPN2"),
  plot_trans_exp_individual(transcript=NULL,class.files$glob_targ_SQ,ExpGenes$targeted_group,"group","EPN2")
)
