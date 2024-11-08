#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for descriptive stats for isoform developmental paper
##
## ---------------------------------

source(paste0(LOGEN,"transcriptome_stats/summarise_classfiles.R"))

## ------- whole transcriptome -------

wholetarg_annotatedGenes_summary_table <- descriptives_summary(class.files$glob_SQ, class.files$glob_SQ_annoGene)
wholecollapsed_annotatedGenes_summary_table <- descriptives_summary(class.files$glob_collapsed, class.files$glob_collapsed_annoGene)

# abundance of novel transcripts vs known transcripts of known genes 
# sum the mean of the counts across all the whole samples
# t-test 
novelMean <- annoGenesStats$novelTrans %>% select(contains("Whole")) %>% apply(., 1, mean) 
knownMean <- annoGenesStats$annoTrans %>% select(contains("Whole")) %>% apply(., 1, mean) 
dat <- rbind(reshape2::melt(novelMean) %>% mutate(associated_transcript = "novel"), reshape2::melt(knownMean) %>% mutate(associated_transcript = "known"))
res <- t.test(value ~ associated_transcript, data = dat)
res$p.value
class.files$glob_SQ_annoGene %>% filter(structural_category %in% c("NIC","NNC")) %>% mutate(FL = preReads + postReads) %>% arrange(-FL)

# proteogenomics pipeline
proteinInput$t2p.collapse.AnnoGenes <- proteinInput$t2p.collapse %>% filter(base_acc %in% class.files$glob_SQ_annoGene$isoform)
message("Number of RNA transcripts: ", length(unique(proteinInput$t2p.collapse.AnnoGenes$pb_accs)))
message("Percentage of RNA transcripts: ", length(unique(proteinInput$t2p.collapse.AnnoGenes$pb_accs))/nrow(class.files$glob_SQ_annoGene))
message("Number of RNA isoforms: ", length(unique(proteinInput$t2p.collapse.AnnoGenes$base_acc)))
Cpat$whole.AnnoGenes <- Cpat$whole %>% filter(seq_ID %in% proteinInput$t2p.collapse.AnnoGenes$base_acc)
message("Mean ORF size: ", mean(Cpat$whole.AnnoGenes$ORF), " (", sd(Cpat$whole.AnnoGenes$ORF),")")
message("Number of coding ORFs", nrow(Cpat$whole.AnnoGenes %>% filter(Coding_prob >=0.364)))
Cpat$whole.AnnoGenes.coding <- Cpat$whole.AnnoGenes %>% filter(Coding_prob >=0.364) %>% filter(seq_ID %in% proteinInput$filtered$pb)
message("Number of coding ORFs with stop codons i.e. protein isoform: ", nrow(Cpat$whole.AnnoGenes.coding))
Cpat$whole.AnnoGenes.noncoding <- Cpat$whole.AnnoGenes %>% filter(Coding_prob < 0.364) %>% filter(seq_ID %in% proteinInput$filtered$pb)

# nmd
noncoding_nmd <- proteinInput$nmd %>% filter(pb %in% Cpat$whole.AnnoGenes.noncoding$seq_ID & is_nmd == "True")
coding_nmd <- proteinInput$nmd %>% filter(pb %in% Cpat$whole.AnnoGenes.coding$seq_ID & is_nmd == "True")
nrow(noncoding_nmd)
nrow(coding_nmd)
binom.test(nrow(noncoding_nmd), nrow(Cpat$whole.AnnoGenes.noncoding))
binom.test(nrow(coding_nmd), nrow(Cpat$whole.AnnoGenes.coding), p = 0.5)


# novel transcripts of known genes
proteinInput$t2p.collapse.novelTranscriptsAnnoGenes <- proteinInput$t2p.collapse %>% filter(base_acc %in% annoGenes_stats$novelTrans$isoform)
message("Number of RNA transcripts: ", length(unique(proteinInput$t2p.collapse.novelTranscriptsAnnoGenes$pb_accs)))
message("Number of RNA isoforms: ", length(unique(proteinInput$t2p.collapse.novelTranscriptsAnnoGenes$base_acc)))

# cpat 
Cpat$whole.novelTranscriptsAnnoGenes <- Cpat$whole %>% filter(seq_ID %in% proteinInput$t2p.collapse.novelTranscriptsAnnoGenes$base_acc)
message("Mean ORF size: ", mean(Cpat$whole.novelTranscriptsAnnoGenes$ORF), " (", sd(Cpat$whole.novelTranscriptsAnnoGenes$ORF),")")
message("Number of coding ORFs i.e. protein isoform: ", nrow(Cpat$whole.novelTranscriptsAnnoGenes %>% filter(Coding_prob >=0.364)))

# raw read counts 
rawCountsWhole <- rawCounts %>% filter(grepl("Whole", sample))
removeWholeSample <- rawCountsWhole[rawCountsWhole$group == "Postnatal" & rawCountsWhole$age < 2,"sample"]
rawCountsWhole <- rawCountsWhole[!rawCountsWhole$sample %in% removeWholeSample,]
t.test(counts ~ group, rawCountsWhole)

# prenatal vs postnatal
prevspostlengths <- bind_rows(
  data.frame(lengths = class.files$glob_SQ_annoGene_prenatal$length, group = "prenatal"),
  data.frame(lengths = class.files$glob_SQ_annoGene_postnatal$length, group = "postnatal"),
)
prevspostexons<- bind_rows(
  data.frame(exons = class.files$glob_SQ_annoGene_prenatal$exons, group = "prenatal"),
  data.frame(exons = class.files$glob_SQ_annoGene_postnatal$exons, group = "postnatal"),
)
ggplot(prevspostlengths , aes(x = group, y = lengths)) + geom_boxplot()
t.test(lengths ~ group, data = prevspostlengths)
t.test(exons ~ group, data = prevspostexons)
class.files$glob_SQ_annoGene_prenatal$length
mean(class.files$glob_SQ_annoGene_prenatal$length)
mean(class.files$glob_SQ_annoGene_postnatal$length)
mean(class.files$glob_SQ_annoGene_prenatal$exons)
mean(class.files$glob_SQ_annoGene_postnatal$exons)

# comparison of number of unique and common transcripts in prenatal vs postnatal
commonDevTranscripts <- intersect(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"],
                                  class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"])
message("Number of transcripts to known genes detected in both prenatal and postnatal: ", length(commonDevTranscripts))
message("Number of transcripts detected in prenatal: ",length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"]))
message("% of all prenatal detected: ", length(commonDevTranscripts)/length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"]) * 100)
message("Number of transcripts detected in postnatal: ",length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"]))
message("% of all postnatal detected: ", length(commonDevTranscripts)/length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"]) * 100)
message("Number of transcripts unique to prenatal:", length(setdiff(class.files$glob_SQ_annoGene_prenatal$isoform, class.files$glob_SQ_annoGene_postnatal$isoform)))
message("Number of transcripts unique to postnatal:", length(setdiff(class.files$glob_SQ_annoGene_postnatal$isoform, class.files$glob_SQ_annoGene_prenatal$isoform)))
message("% of unique postnatal to all:", length(setdiff(class.files$glob_SQ_annoGene_postnatal$isoform, class.files$glob_SQ_annoGene_prenatal$isoform))/nrow(class.files$glob_SQ_annoGene) * 100)
message("% of unique prenatal to all:",length(setdiff(class.files$glob_SQ_annoGene_prenatal$isoform, class.files$glob_SQ_annoGene_postnatal$isoform))/nrow(class.files$glob_SQ_annoGene) * 100)

# most abundant transcript in prenatal and not in postnatal
class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads == 0,] %>% arrange(-preReads) %>% .[1,]
class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads == 0,] %>% arrange(-postReads) %>% .[1,]

message("Number of differentially expressed genes: ", length(unique(WholeDESeqGeneSig$age$associated_gene)))

XEscape <- merge(WholeDTE$sex, XExcapeList, by.x = "associated_gene", by.y = "GeneName")
plot_volcano(diff_results=WholeDESeq2$age,stats = TRUE)
nrow(WholeDESeq2Sig$age %>% group_by(associated_gene) %>% tally())
WholeDTE$age %>% group_by(structural_category) %>% tally()
17308/85428
63204/85428

## ------- long-read sequencing datasets comparisons -------

# number of overlaps between Leung et al. 2021 and whole dataset

comparison_dataset <- function(gfftmap, classfiles, altDataset){

        Detected_alt <- unique(gfftmap[gfftmap$class_code == "=","qry_id"])
        Detected_RB <- unique(gfftmap[gfftmap$qry_id %in% Detected_alt,"ref_id"])
        Unique_alt <- setdiff(unique(gfftmap[gfftmap$class_code != "=","qry_id"]),unique(gfftmap[gfftmap$class_code == "=","qry_id"]))
        
        # subset to the classfiles
        Unique_RB <- classfiles[!classfiles$isoform %in% Detected_RB,]
        Common_RB <- classfiles[classfiles$isoform %in% Detected_RB,]

        message("Total number of transcripts alternative long read sequencing dataset: ",length(altDataset$isoform))
        message("Number overlap from alternative long read sequencing dataset: ",length(Detected_alt))
        message("Percentage covered from alternative long read sequencing dataset: ", round(length(Detected_alt)/length(unique(altDataset$isoform)) * 100,2))
        message("Number overlap from Bamford dataset: ", length(unique(Common_RB$isoform)))
        message("Percentage overlap from Bamford dataset: ", round(length(unique(Common_RB$isoform))/nrow(classfiles)* 100,2))
        message("Number unique to Bamford dataset: ", nrow(Unique_RB))

}

comparison_dataset(gfftmapComparisons$cellReports, class.files$glob_SQ, humanCTX)
comparison_dataset(gfftmapComparisons$directRNA, class.files$glob_SQ, directRNA)
comparison_dataset(gfftmapComparisons$BDRNatureComms, class.files$glob_targ_SQ_20AD, BDRNatureComms)
comparison_dataset(gfftmapComparisons$Patowary, class.files$glob_SQ, PatowaryCTX)

ExpressionDiffPacBio <- rbind(Unique_RB %>% mutate(dataset = "Unique") %>% select(nreads, dataset),
                              Common_RB %>% mutate(dataset = "Common") %>% select(nreads, dataset))
res <- wilcox.test(nreads ~ dataset, ExpressionDiffPacBio, exact = FALSE)
res
res$p.value
format(res$p.value, scientific = TRUE)



## ------- differential transcript expression -------

# development
nrow(WholeDTE$age) == length(unique(WholeDTE$age$isoform))
message("Number of differentially expressed isoforms: ", nrow(WholeDTE$age))
message("Number of genes with differentially expressed isoforms: ", length(unique(WholeDTE$age$associated_gene)))
message("Number of upregulated transcripts in postnatal vs prenatal: ", nrow(WholeDTE$age[WholeDTE$age$dirAcrossDev == "upregulated",]))
message("top-ranked transcript differentially expressed across development")
head(WholeDTE$age %>% arrange(padj))

# binomial test of enrichment of upregulated transcripts in post-natal
binom <- binom.test(nrow(WholeDTE$age[WholeDTE$age$dirAcrossDev == "upregulated",]), nrow(WholeDTE$age))
print(binom$p.value)

# binomial test of novel DETs
WholeDTE$ageNovelTranscripts <- WholeDTE$age %>% filter(grepl("novel", associated_transcript))
message("Number of differentially expressed novel isoforms: ", nrow(WholeDTE$ageNovelTranscripts), "(", nrow(WholeDTE$ageNovelTranscripts)/nrow(WholeDTE$age),")")
binom <- binom.test(nrow(WholeDTE$ageNovelTranscripts), nrow(WholeDTE$age))
print(binom$p.value)
message("top-ranked novel transcript differentially expressed across development")
head(WholeDTE$age %>% arrange(padj) %>% filter(associated_transcript == "novel"))

# Gene ontology input for genes with top-ranked DETs
View(unique(WholeDTE$age$associated_gene)[1:100])

# antisense DETs
message("Number of antisense DETs: ", nrow(WholeDTE$age[WholeDTE$age$structural_category == "Antisense",]))
message("% of antisense DETs: ", nrow(WholeDTE$age[WholeDTE$age$structural_category == "Antisense",])/nrow(WholeDTE$age) * 100)
View(WholeDTE$age[WholeDTE$age$structural_category == "Antisense" & WholeDTE$age$exons > 1,])

## sex
message("Number of differentially expressed transcrips by sex: ", nrow(WholeDTE$sex))
message("Number of genes with differentially expressed transcrips by sex: ", length(unique(WholeDTE$sex$associated_gene)))
message("Number of transcripts on the X and Y chromosome: ", nrow(WholeDTE$sex %>% filter(chrom %in% c("chrX","chrY"))))
head(WholeDTE$sex %>% filter(!chrom %in% c("chrX","chrY")) %>% arrange(padj))
message("common transcripts differentially expressed across development and by sex")
intersect(WholeDTE$sex$isoform, WholeDTE$age$isoform)
head(WholeDTE$sex %>% filter(chrom %in% c("chrX","chrY")) %>% arrange(padj))

# targeted transcriptome
message("Mean number of exons, (sd): ", round(mean(class.files$targ_SQ$exons),2)," (", round(sd(class.files$targ_SQ$exons),2),")")
length(unique(class.files$targ_SQ$associated_gene))
nrow(class.files$targ_SQ)
mean(class.files$targ_SQ$length)
sd(class.files$targ_SQ$length)
min(class.files$targ_SQ$length)
max(class.files$targ_SQ$length)
numIsogeneTallyTargeted = class.files$targ_SQ%>% group_by(associated_gene) %>% tally()
min(numIsogeneTallyTargeted$n)
max(numIsogeneTallyTargeted$n)
median(numIsogeneTallyTargeted$n)
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 1,])
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 1,])/length(unique(numIsogeneTallyTargeted$associated_gene))
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 10,])
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 10,])/length(unique(numIsogeneTallyTargeted$associated_gene))
nrow(class.files$targ_SQ[class.files$targ_SQ$structural_category == "NNC",])/nrow(class.files$targ_SQ)

plot_volcano(diff_results=TargetedDESeq2$age,stats = TRUE)
nrow(TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n>10))
TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% arrange(-n)
159/length(unique(TargetedDESeq2Sig$age$associated_gene))
TargetedDESeq2Sig$age %>% group_by(structural_category) %>% tally()

TargetedDESeq2Sig$age %>% filter(associated_gene %in% c(TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n==1) %>% .$associated_gene))
TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n==1)
807/nrow(TargetedDESeq2Sig$age)
4692/nrow(TargetedDESeq2Sig$age)


# list of genes escaping X inactivation in DTE by sex
WholeDTE$sex
XExcape$GeneName


## ---------- Alternative splicing events ----------

FICLE_class_final_exp <- FICLE_class_final_exp %>% mutate(A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3),
                                                          ES = as.numeric(ES), ESpre = ES * normPre/sum(ES), ESpost = ES * normPost/sum(ES),
                                                          IR = as.numeric(IR), IRpre = IR * normPre/sum(IR), IRpost = IR * normPost/sum(IR),
                                                          A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3))

FICLEES <- FICLE_class_final_exp %>% select(isoform, ESpre, ESpost) %>% reshape2::melt(variable.name = "development", value.name = "normES") 
FICLEIR <- FICLE_class_final_exp %>% select(isoform, IRpre, IRpost) %>% reshape2::melt(variable.name = "development", value.name = "normIR") 
FICLEA5A3 <- FICLE_class_final_exp %>% select(isoform, A5A3pre, A5A3post) %>% reshape2::melt(variable.name = "development", value.name = "normA5A3") 

ggplot(FICLEES, aes(x = development, y = log10(normES))) + geom_boxplot()
ggplot(FICLEIR, aes(x = development, y = log10(normIR))) + geom_boxplot()
ggplot(FICLEA5A3, aes(x = development, y = log10(normA5A3))) + geom_boxplot()


t.test(normES ~ development, FICLEES)
t.test(normA5A3 ~ development, FICLEA5A3)
t.test(normIR ~ development, FICLEIR)

row.names(FICLE_class) <- FICLE_class$isoform
FICLE_class <- FICLE_class %>% mutate_if(is.character, as.numeric)
FICLE_class <- FICLE_class[complete.cases(FICLE_class), ]


count = 1
AS_events <- data.frame()
for(i in colnames(FICLE_class)){
  AS_events[count,1] <- i
  AS_events[count,2] <- sum(FICLE_class[[i]])
  count <- count + 1
}
colnames(AS_events) <- c("AS_event","Number")
ggplot(AS_events, aes(x = reorder(AS_event,-Number), y = Number)) + geom_bar(stat = "identity")



## ---------- Proteogenomics ----------

message("Number of RNA Transcripts to all genes:", nrow(class.files$glob_SQ))
message("Number of coding RNA isoforms:", length(unique(proteinInput$t2p.collapse.refined$pb_acc)))
message("Number of coding RNA isoforms:", length(unique(proteinInput$t2p.collapse.refined$corrected_acc)))
length(unique(class.files$ptarg_filtered$corrected_acc)) # -1 to not include "NA" 

setdiff(class.files$glob_SQ$isoform, c(as.character(proteinInput$no_orf$V1),as.character(proteinInput$cpat$seq_ID)))
# these are the isoforms that have been filtered from monoexonic isoforms from multiexonic genes
setdiff(c(as.character(proteinInput$no_orf$V1),as.character(proteinInput$cpat$seq_ID)),class.files$glob_SQ$isoform)

message("Number of transcripts with noORF: ", length(as.character(proteinInput$no_orf$V1)))
message("Number of transcripts with ORF: ", length(as.character(proteinInput$cpat$seq_ID)))
Rank1 <- proteinInput$mapped[proteinInput$mapped$orf_rank == "1",]
message("Number of transcripts with ORF ranked 1 and has stop codons: ", length(unique(Rank1[Rank1$has_stop_codon == "TRUE","transcript_id"])))
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0: ", length(unique(proteinInput$t2p.collapse.refined$pb_accs)))

RefinedORF <- proteinInput$cpat[proteinInput$cpat$seq_ID %in% proteinInput$t2p.collapse.refined$pb_accs,]
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0.364: ", nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,]))
nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,])/nrow(class.files$glob_SQ)

RefinedORFCoding <- RefinedORF[RefinedORF$Coding_prob > 0.364,]
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0.364, collapsed: ",
        length(unique(proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"corrected_acc"])))
590504/nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,])

GSselectedcollapsedID <- proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"base_acc"]
CorrectedcollapsedID <- proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"corrected_acc"]

# note one transcript was filtered from SQANTI protein so off by 1 when calculating difference
message("Number of transcripts with ORF...collapsed: ", length(unique(setdiff(RefinedORFCoding$seq_ID, CorrectedcollapsedID))))
ORFCorrectedCollapsedID <- unique(setdiff(RefinedORFCoding$seq_ID, CorrectedcollapsedID))
class.files$targ_filtered %>% filter(isoform %in% ORFCorrectedCollapsedID) %>% group_by(structural_category) %>% tally()

class.files$protein_filtered_final <- class.files$protein_filtered[class.files$protein_filtered$pb %in% GSselectedcollapsedID,]
message("Number of protein products after filtering: ", nrow(class.files$protein_filtered_final))

## ---------- Mass spectrometry ----------

# All peptides
peptidesFiltered <- lapply(peptides, function(x) x %>% dplyr::select('Base Sequence', 'Full Sequence', 'Protein Accession'))
peptidesFiltered <- dplyr::bind_rows(peptidesFiltered)

# remove duplicated rows (due to merging samples)
peptidesFiltered <- peptidesFiltered %>% distinct()
message("Total number of peptides:", length(unique(peptidesFiltered$`Base Sequence`)))
#do not use nrows for peptide counts #peptidesFiltered[peptidesFiltered$`Base Sequence` == "AADAEAEVASLNR",]

# determine number of transcripts validated in terms of splice junction (split protein accession column)
accession_split <- lapply(peptidesFiltered[["Protein Accession"]], function(x) strsplit(x, "|", fixed = TRUE)[[1]])
accession_split <- unique(unlist(accession_split))

# remove decoy protein accession
accession_split <- accession_split[!grepl("DECOY", accession_split)]
accession_split <- accession_split[grepl("ONT", accession_split)]
message("Total number of trancripts with mass-spec validation:", nrow(class.files$glob_SQ[class.files$glob_SQ$isoform %in% accession_split,]))
message("Total number of novel trancripts with mass-spec validation:", nrow(class.files$glob_SQ[class.files$glob_SQ$isoform %in% accession_split,] %>% filter(!structural_category %in% c("FSM", "ISM"))))

## novel peptides
novelPeptides <- bind_rows(novelPeptides)
novelPeptides <- novelPeptides[novelPeptides$acc %in% class.files$glob_SQ$isoform,] %>% select(gene, acc, seq)
novelPeptides <- novelPeptides %>% distinct()
#novelPeptides[novelPeptides$acc == "ONT12_1620_26499",]

novelPeptides <- novelPeptides %>% mutate(length = nchar(as.character(seq)))
novelPeptidesunique <- novelPeptides %>% group_by(acc) %>% top_n(1, length)
colnames(novelPeptides) <- c("gene","isoform","protein_seq")
novelPeptides <- merge(novelPeptides,class.files$glob_SQ, by = "isoform") %>% select(isoform, associated_gene, structural_category, protein_seq)
novelPeptidesTranscripts <- novelPeptides[!novelPeptides$structural_category %in% c("FSM","ISM"),]
write.table(novelPeptidesTranscripts, paste0(dirnames$output,"novelTranscriptsPeptides.txt"), sep = "\t", quote = F)


## ---------- bambu ----------

message("Number of total transcripts: ", nrow(class.files$bambu))
message("Number of all genes: ", length(unique((class.files$bambu$associated_gene))))
message("Number of annotated known genes: ", length(unique((class.files$bambu_annoGene$associated_gene))))
message("Number of total transcripts to known genes: ", nrow(class.files$bambu_annoGene))
message("Number of novel transcripts to annotated known genes: ", nrow(bambuAnnoGeneStats$novelTrans), "( ", 
        round(nrow(bambuAnnoGeneStats$novelTrans)/nrow(class.files$bambu_annoGene) * 100,2), "%)")

message("Mean length (sd): ", round(mean(class.files$bambu_annoGene$length),2)," (", round(sd(class.files$bambu_annoGene$length),2),")")
message("Min - max length: ", round(min(class.files$bambu_annoGene$length),2)," - ", round(max(class.files$bambu_annoGene$length),2),"")
message("Mean number of exons, (sd): ", round(mean(class.files$bambu_annoGene$exons),2)," (", round(sd(class.files$bambu_annoGene$exons),2),")")


# most isomorphic gene
message("Number of transcripts detected to HNRNPK: ", nrow(class.files$bambu_annoGene[class.files$bambu_annoGene$associated_gene == "HNRNPK",]))

# most isoformic novel transcript
class.files$glob_targ_SQ_counts[class.files$glob_targ_SQ_counts$associated_transcript == "novel",] %>% arrange(-nreads) %>% select(isoform, structural_category, exons, associated_gene, associated_transcript, nreads, nsamples)
class.files$bambu[class.files$bambu$associated_gene == "RPS27A",]

novelPeptidesTranscripts <- read.table(paste0(dirnames$output,"novelTranscriptsPeptides.txt"))
class.files$glob_targ_SQ_counts[class.files$glob_targ_SQ_counts$isoform %in% novelPeptidesTranscripts$isoform,] %>% arrange(-nreads) 
class.files$bambu[class.files$bambu$associated_gene == "HMOX2",]
  
