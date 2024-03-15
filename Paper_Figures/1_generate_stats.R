#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for descriptive stats for isoform developmental paper
##
## ---------------------------------

## ------- whole transcriptome -------

# number of transcripts and genes to annotated genes
message("Number of total transcripts to all genes: ", nrow(class.files$glob_SQ_annoGene))
message("Number of annotated known genes: ", length(unique((class.files$glob_SQ_annoGene$associated_gene))))

# length summary
message("Mean length (sd): ", round(mean(class.files$glob_SQ_annoGene$length),2)," (", round(sd(class.files$glob_SQ_annoGene$length),2),")")
message("Min - max length: ", round(min(class.files$glob_SQ_annoGene$length),2)," - ", round(max(class.files$glob_SQ_annoGene$length),2),"")

# number of transcripts summary
numIsogene = numIsoGene(class.files$glob_SQ_annoGene,stats=TRUE)
numIsogeneTally = class.files$glob_SQ_annoGene %>% group_by(associated_gene) %>% tally()
message("Most isomorphic gene: ")
#numIsogeneTally %>% arrange(-n)
message("Mean number of isoforms (sd): ", round(mean(numIsogeneTally$n),2)," (", round(sd(numIsogeneTally$n),2),")")
message("Min - max number of isoforms: ", round(min(numIsogeneTally$n),2)," - ", round(max(numIsogeneTally$n),2))
nrow(numIsogeneTally[numIsogeneTally$n >= 10,])
nrow(numIsogeneTally[numIsogeneTally$n >= 10,])/length(unique(numIsogeneTally$associated_gene))

# exon summary 
message("Mean number of exons, (sd): ", round(mean(class.files$glob_SQ_annoGene$exons),2)," (", round(sd(class.files$glob_SQ_annoGene$exons),2),")")
meanExonGene <- aggregate(class.files$glob_SQ_annoGene[,"exons"], list(class.files$glob_SQ_annoGene$associated_gene), mean)
message("Mean number of exons for any given gene (sd): ", round(mean(meanExonGene$x),2)," (", round(sd(meanExonGene$x),2),")")

# number of novel and known transcripts
message("Number of total transcripts to annotated known genes: ", nrow(class.files$glob_SQ_annoGene))

annoGenesStats <- list(
  novelTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel",],
  annoTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript != "novel",],
  NIC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NIC",],
  NNC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NNC",]
)

message("Number of novel transcripts to annotated known genes: ", nrow(annoGenesStats$novelTrans), "( ", 
        round(nrow(annoGenesStats$novelTrans)/nrow(class.files$glob_SQ_annoGene) * 100,2), "%)")

message("Number of annotated known genes with novel transcripts: ", length(unique(annoGenesStats$novelTrans$associated_gene)), "( ", 
        round(length(unique(annoGenesStats$novelTrans$associated_gene))/length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100,2), "%)")

message("Mean length (sd): ", round(mean(annoGenesStats$novelTrans$length),2)," (", round(sd(annoGenesStats$novelTrans$length),2),")")
message("Mean number of exons (sd): ", round(mean(annoGenesStats$novelTrans$exons),2)," (", round(sd(annoGenesStats$novelTrans$exons),2),")")


message("Number of known transcripts to annotated known genes: ", nrow(annoGenesStats$annoTrans), "( ", 
        round(nrow(annoGenesStats$annoTrans)/nrow(class.files$glob_SQ_annoGene) * 100,2), "%)")

message("Number of NIC transcripts to annotated known genes: ", nrow(annoGenesStats$NIC), "( ", 
        round(nrow(annoGenesStats$NIC)/nrow(annoGenes_stats$novelTrans) * 100,2), "%)")

message("Number of NNC transcripts to annotated known genes: ", nrow(annoGenesStats$NNC), "( ", 
        round(nrow(annoGenesStats$NNC)/nrow(annoGenes_stats$novelTrans) * 100,2), "%)")

# abundance of novel transcripts vs known transcripts of known genes 
# sum the mean of the counts across all the whole samples
# t-test 
novelMean <- annoGenesStats$novelTrans %>% select(contains("Whole")) %>% apply(., 1, mean) 
knownMean <- annoGenesStats$annoTrans %>% select(contains("Whole")) %>% apply(., 1, mean) 
dat <- rbind(reshape2::melt(novelMean) %>% mutate(associated_transcript = "novel"), reshape2::melt(knownMean) %>% mutate(associated_transcript = "known"))
t.test(value ~ associated_transcript, data = dat)
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
message("prenatal vs postnatal")
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


commonDevTranscripts <- intersect(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"],class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"])
message("Number of transcripts to known genes detected in both prenatal and postnatal: ", length(commonDevTranscripts))
message("% of all prenatal detected: ", length(commonDevTranscripts)/length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,"isoform"]) * 100)
message("% of all postnatal detected: ", length(commonDevTranscripts)/length(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,"isoform"]) * 100)

message("Number of transcripts unique to postnatal:", length(setdiff(class.files$glob_SQ_annoGene_postnatal$isoform, class.files$glob_SQ_annoGene_prenatal$isoform)))
message("Number of transcripts unique to prenatal:", length(setdiff(class.files$glob_SQ_annoGene_prenatal$isoform, class.files$glob_SQ_annoGene_postnatal$isoform)))

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


# number of overlaps between Leung et al. 2021 and whole dataset
message("Number overlap: ",length(unique(gfftmap[gfftmap$class_code == "=","qry_id"])))
message("Percentage overlap: ",length(unique(gfftmap[gfftmap$class_code == "=","qry_id"]))/nrow(humanCTX))


## ------- differential transcript expression -------

# development
message("Number of differentially expressed isoforms: ", nrow(WholeDTE$age))
message("Number of genes with differentially expressed isoforms: ", length(unique(WholeDTE$age$associated_gene)))
message("Number of upregulated transcripts in postnatal vs prenatal: ", nrow(WholeDTE$age[WholeDTE$age$dirAcrossDev == "upregulated",]))
table(WholeDTE$age$structural_category)

# binomial test of enrichment of upregulated transcripts in post-natal
binom <- binom.test(nrow(WholeDTE$age[WholeDTE$age$dirAcrossDev == "upregulated",]), nrow(WholeDTE$age))
print(binom$p.value)

# binomial test of novel DETs
WholeDTE$ageNovelTranscripts <- WholeDTE$age %>% filter(grepl("novel", associated_transcript))
message("Number of differentially expressed novel isoforms: ", nrow(WholeDTE$ageNovelTranscripts), "(", nrow(WholeDTE$ageNovelTranscripts)/nrow(WholeDTE$age),")")
binom <- binom.test(nrow(WholeDTE$ageNovelTranscripts), nrow(WholeDTE$age))
print(binom$p.value)

# Gene ontology input for genes with top-ranked DETs
View(unique(WholeDTE$age$associated_gene)[1:100])

# antisense DETs
message("Number of antisense DETs: ", nrow(WholeDTE$age[WholeDTE$age$structural_category == "Antisense",]))
View(WholeDTE$age[WholeDTE$age$structural_category == "Antisense" & WholeDTE$age$exons > 1,])

## sex
message("Number of differentially expressed transcrips by sex: ", nrow(WholeDTE$sex))
message("Number of genes with differentially expressed transcrips by sex: ", length(unique(WholeDTE$sex$associated_gene)))

nrow(WholeDTE$sex %>% filter(!chrom %in% c("chrX","chrY")))
nrow(unique(WholeDTE$sex %>% filter(!chrom %in% c("chrX","chrY")) %>% .["associated_gene"]))

intersect(WholeDTE$sex$isoform, WholeDTE$age$isoform)

class.files$targ_SQ %>% group_by(associated_gene) %>% tally() %>% arrange(-n)

Cpat$whole

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
