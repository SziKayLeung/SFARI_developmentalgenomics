#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for supplementary figures
##
## ---------------------------------


SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs")

# age
pAge <- ages(phenotype$WholeTargeted)

# sensitivity curve of filtering in targeted dataset
no_of_isoforms_sample(class.files$targ_SQ)

# comparison of whole vs targeted datasets across matched samples
comp = whole_vs_targeted_plots(classfiles=class.files$glob_targ_SQ_counts, wholeSamples=wholematchedsamples, targetedSamples=targetedmatchedsamples, targetGene=selectedTargetGenes)
pdf(paste0(output_dir,"/UniqueIsoformsWholeDataset.pdf"), width = 10, height = 30)
comp[[1]]
plot_grid(comp[[3]])
dev.off()

plot_grid(plotlist = comp[1:6], labels = c("A","B","C","D","E","F"))

plot_cupcake_collapse_sensitivity(class.files$targ_SQ,"All target genes")


# Number of DTEs in development by structural category
WholeDTE$age %>% group_by(structural_category, dirAcrossDev) %>% tally() %>% 
  filter(structural_category != "NA") %>%
  ggplot(., aes(x = reorder(structural_category,n), y = n, fill = dirAcrossDev)) + geom_bar(stat = "identity",position = position_dodge()) + 
  coord_flip() +
  theme_classic() + labs(y = "Number of transcripts", x = "Structural category") +
  scale_fill_discrete(name = "Direction across development") + theme(legend.position = "top")

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

## --- Comparison between PacBio and fetal dataset 

PacBioWhole <- rbind(
  humanCTX[humanCTX$isoform %in% DetectedBoth,c("structural_category","totalFL","exons","length")] %>% mutate(Dataset = "Detected"),
  humanCTX[humanCTX$isoform %in% Unique,c("structural_category","totalFL","exons","length")] %>% mutate(Dataset = "Unique")
) 
pPacBioWhole <- list(
  ggplot(PacBioWhole, aes(x = Dataset, y = log10(totalFL))) + geom_boxplot() + theme_classic() + labs(x = "", y = "log10 FL read count"),
  ggplot(PacBioWhole, aes(x = Dataset, y = exons)) + geom_boxplot() + theme_classic() + labs(x = "", y = "Number of exons"),
  ggplot(PacBioWhole, aes(x = Dataset, y = length)) + geom_boxplot() + theme_classic() + labs(x = "", y = "Length (bp)")
)

plot_grid(plotlist = pPacBioWhole, labels = c("A","B","C"))

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

tallyReads <- class.files$glob_targ_SQ_counts %>% filter(!grepl("novel", associated_gene)) %>% 
  group_by(nreads, nsamples) %>% 
  tally()
tallyReads <- as.data.frame(tallyReads) %>% mutate(perc = n/sum(n) * 100)


tallyReads <- tallyReads %>%
  arrange(nreads, nsamples) %>%  # Sort the data (optional, depending on the desired order)
  mutate(cum_sum = cumsum(n),
         cum_percentage = cum_sum / sum(n) * 100)

tallyReads <- tallyReads %>% mutate(label = ifelse(cum_percentage < 70, paste0(nreads, "reads,", nsamples, "samples"), ""))
ggplot(tallyReads, aes(x = nreads, y = cum_percentage, colour = factor(nsamples), label = label)) + geom_point() +
  scale_x_continuous(breaks = seq(min(tallyReads$nreads), max(tallyReads$nreads), by = 50000)) +
  geom_label_repel(show.legend = FALSE) +
  theme_classic() +
  labs(x = "Number of reads", y = "Cumulative percentage of transcripts annotated to known genes", colour = "Number of samples")


## ----- correlation of RIN with number of transcripts -------

datWholeCounts <- class.files$glob_targ_SQ %>% select(contains("Whole"),-"whole_nsamples",-"whole_nreads")
## to check below command about colSums
#length(datWholeCounts[,"Whole11831"][datWholeCounts[,"Whole11831"] != 0])
#length(datWholeCounts[,"Whole11831"][datWholeCounts[,"Whole11831"] == 0])
#nrow(datWholeCounts) == length(datWholeCounts[,"Whole11831"][datWholeCounts[,"Whole11831"] != 0]) + length(datWholeCounts[,"Whole11831"][datWholeCounts[,"Whole11831"] == 0])

# tally the occurences in each column, not equal to 0 i.e. a read detected for transcript, and therefore transcript detected in dataset
TranscriptPerSample <- colSums(datWholeCounts  != 0) %>%  reshape2::melt(., value.name = "counts") %>% tibble::rownames_to_column(., var = "ID")

# merge with manifest to get the RIN numbers
TranscriptPerSample <- merge(TranscriptPerSample, manifest, by = "ID")

# postnatal dataset
PostnatalTranscriptPerSample <- TranscriptPerSample[TranscriptPerSample$Group == "Postnatal",]
cor.test(PostnatalTranscriptPerSample$RIN, PostnatalTranscriptPerSample$counts)

# postnatal dataset (with RIN > 5)
PostnatalTranscriptPerSampleGoodRIN <- TranscriptPerSample[TranscriptPerSample$Group == "Postnatal" & TranscriptPerSample$RIN >= 4.5,]
cor.test(PostnatalTranscriptPerSampleGoodRIN$RIN, PostnatalTranscriptPerSampleGoodRIN$counts)

cor.test(TranscriptPerSample$RIN, TranscriptPerSample$counts)
pRINTranscripts <- ggplot(TranscriptPerSample, aes(x = RIN, y = counts, colour = Group)) + 
  geom_point(size = 3) +
  labs(x = "RIN", y = "Number of Transcripts") + 
  theme_classic()

