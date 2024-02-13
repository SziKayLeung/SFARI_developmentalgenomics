#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary figures
##
## ---------------------------------


SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs")

# age
pAge <- ages(phenotype$WholeTargeted)

# sensitivity curve of filtering in targeted dataset
no_of_isoforms_sample(class.files$targ_SQ)

# comparison of whole vs targeted datasets across matched samples
comp = whole_vs_targeted_plots(classfiles=class.files$glob_targ_SQ, wholeSamples=wholematchedsamples, targetedSamples=targetedmatchedsamples, targetGene=TargetGene)
pdf(paste0(output_dir,"/UniqueIsoformsWholeDataset.pdf"), width = 10, height = 30)
comp[[7]]
dev.off()

plot_grid(plotlist = comp[1:6], labels = c("A","B","C","D","E","F"))

plot_cupcake_collapse_sensitivity(class.files$targ_SQ,"All target genes")

# prenatal vs postnatal 
ggplot(geneNum, aes(x = prenatalTranscripts, y = postnatalTranscripts)) + geom_point() +
  theme_classic() + labs(x = "Number of transcripts: prenatal", y = "Number of transcripts: postnatal")


# sensitivity plots
#pSensitivity(class.files$targ_SQ)
#ggsave(file=paste0(output_dir,"/cumulativeSensitivityTargeted.png"), dpi=400, width = 20, height = 20, units = "cm")

#pSensitivity(class.files$glob_SQ)
#ggsave(file=paste0(output_dir,"/cumulativeSensitivityWhole.png"), dpi=400, width = 20, height = 20, units = "cm")

## number of prenatal vs postnatal 
## prenatal vs postnatal
devVennNums <- twovenndiagrams(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"],
                               class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"],"Pre-natal","Post-natal")
grid.draw(devVennNums)

class.files$glob_SQ$DevStatus <- apply(class.files$glob_SQ, 1, function(x) identify_dataset_by_counts(x[["postReads"]], x[["preReads"]], "postnatal","prenatal"))
class.files$glob_SQ_annoGene %>% group_by(structural_category,DevStatus) %>% tally(name = "num") %>% 
  ggplot(., aes(x = structural_category, y = num, fill = DevStatus)) + geom_bar(stat = "identity", position="dodge") +
  labs(x = "Structural category", y = "Number of transcripts of annotated genes") +
  scale_fill_discrete(label = c("Both","Post-natal only", "Pre-natal only"), name = NULL) +
  theme_classic() + theme(legend.position = "top")

class.files$glob_SQ_annoGene %>% filter(DevStatus == "prenatal") %>% arrange(-preReads)
class.files$glob_SQ_annoGene %>% filter(DevStatus == "postnatal") %>% arrange(-postReads)
uniquePrenatalGenes <- unique(class.files$glob_SQ_annoGene %>% filter(DevStatus == "prenatal") %>% select(associated_gene))
write.table(uniquePrenatalGenes, paste0(dirnames$output,"uniquePrenatalGenes.txt"), row.names = F, col.names = F, quote = F)


pIF <- list(
  ontNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_ont$normAll,
                                                      pheno=phenotype$targeted_rTg4510_ont,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$ontDIUGeno$keptIso))),
  
  isoNorm = lapply(Targeted$Genes, function(x) plotIF(x,
                                                      ExpInput=Exp$targ_iso$normAll,
                                                      pheno=phenotype$targeted_rTg4510_iso,
                                                      cfiles=class.files$targ_all,
                                                      design="time_series",
                                                      majorIso=row.names(TargetedDIU$isoDIUGeno$keptIso)))
)
for(i in 1:length(pIF)){names(pIF[[i]]) <- Targeted$Genes}


plot_grid(plotlist = pAge, labels = c("A","B","C","D"))


# disease focus
num_disease_focus_DTE(TargetedDESeq2Sig$age,disease_list$SCHEMA$Gene,"SCHEMA")
TargetedDESeq2Sig$age %>% filter(associated_gene %in% disease_list$SCHEMA$Gene) %>% arrange(padj)


## GRIA3
pTargetedDESeq2SigAgeSchema <- list()
count=1
for(t in c("ONTX_7115_8288","ONTX_7115_7279","ONTX_7115_8547")){
  pTargetedDESeq2SigAgeSchema[[count]] <- 
    plot_grid(plot_trans_exp_individual(t,class.files$glob_targ_SQ,Exp$targeted_group,"group"),
              plot_trans_exp_lifetime(t,class.files$glob_targ_SQ,Exp$targeted_group))
  count = count + 1
}

RefIsoforms <- lapply(c("GRIN2A","GRIA3"), function(x) unique(gtf$ref[gtf$ref$gene_name == x & !is.na(gtf$ref$transcript_id), "transcript_id"]))
names(RefIsoforms ) <- c("GRIN2A","GRIA3")
pTargetedDESeq2SigAgeSchemaTracks <- ggTranPlots(gtf$merged,class.files$glob_targ_SQ,
            isoList = c("ONTX_7115_8288","ONTX_7115_7279","ONTX_7115_8547",RefIsoforms$GRIA3),
            colours = c(wes_palette("Cavalcanti1")[5],wes_palette("Darjeeling2")[2],wes_palette("Royal1")[4],rep("#0C0C78",length(RefIsoforms$GRIA3)+1)), 
            lines =  c(wes_palette("Cavalcanti1")[5],wes_palette("Darjeeling2")[2],wes_palette("Royal1")[4],rep("#0C0C78",length(RefIsoforms$GRIA3)+1)), 
            gene = "GRIN2A",simple=TRUE)

pdf(paste0(output_dir,"/targeted_GRIA3.pdf"), width = 15, height = 10)
plot_grid(plot_grid(plotlist = pTargetedDESeq2SigAgeSchema,nrow=3),pTargetedDESeq2SigAgeSchemaTracks,nrow=1,rel_widths = c(0.6,0.4))
dev.off()

pdf(paste0(output_dir,"/targeted_schema.pdf"), width = 10, height = 10)
num_disease_focus_DTE(TargetedDESeq2Sig$age,disease_list$SCHEMA$Gene,"SCHEMA")
dev.off()