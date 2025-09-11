#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for descriptive stats for isoform developmental paper
##
## ---------------------------------

source(paste0(LOGEN,"/transcriptome_stats/summarise_classfiles.R"))

## ------- whole transcriptome -------

wholetarg_annotatedGenes_summary_table <- descriptives_summary(class.files$glob_SQ, class.files$glob_SQ_annoGene)
wholecollapsed_annotatedGenes_summary_table <- descriptives_summary(class.files$glob_collapsed, class.files$glob_collapsed_annoGene)
wholecollapsed_annotatedGenes_sqanti_default_summary_table <- descriptives_summary(class.files$glob_SQ_default, class.files$glob_SQ_default)

View(class.files$glob_SQ_default %>% group_by(structural_category) %>% tally())
View(class.files$glob_SQ %>% group_by(structural_category) %>% tally())
View(class.files$glob_SQ %>% group_by(structural_category) %>% tally())


# HNRPNK
class.files$glob_SQ %>% group_by(associated_gene) %>% tally() %>% arrange(-n) %>% head(.)
message("Number of isoforms annotated to HNRPNK: ", nrow(class.files$glob_SQ_annoGene %>% filter(associated_gene == "HNRNPK")))

# abundance of novel transcripts vs known transcripts of known genes 
# sum the mean of the counts across all the whole samples
# t-test 

normalized_glob_SQ <- class.files$glob_SQ %>% tibble::column_to_rownames(., var = "isoform") %>%
  select(contains("Whole", ignore.case = FALSE)) %>%   # Select columns containing "Whole"
  mutate(across(everything(), ~ .x / sum(.x, na.rm = TRUE) * 1000000))  # Normalize each column

annotatedGenesNovelTranscripts <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel","isoform"]
annotatedGenesKnownTranscripts <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript != "novel","isoform"]
nrow(annoGenesStats$novelTrans) == length(annotatedGenesNovelTranscripts)
nrow(annoGenesStats$annoTrans) == length(annotatedGenesKnownTranscripts)

novelMean <- normalized_glob_SQ[annotatedGenesNovelTranscripts,] %>% apply(., 1, mean) 
knownMean <- normalized_glob_SQ[annotatedGenesKnownTranscripts,] %>% apply(., 1, mean) 
dat <- rbind(reshape2::melt(novelMean) %>% mutate(associated_transcript = "novel"), reshape2::melt(knownMean) %>% mutate(associated_transcript = "known"))
dat %>% group_by(associated_transcript) %>% summarise(median = mean(value))
res <- t.test(value ~ associated_transcript, data = dat)
res$p.value

message("Percentage of novel transcripts with less than 5 reads:", nrow(annoGenesStats$novelTrans %>% filter(whole_nreads <= 5))/nrow(annoGenesStats$novelTrans))
message("Percentage of novel transcripts with less than 5 reads:", nrow(annoGenesStats$annoTrans %>% filter(whole_nreads <= 5))/nrow(annoGenesStats$annoTrans))
class.files$glob_SQ_annoGene %>% filter(structural_category %in% c("NIC","NNC")) %>% mutate(FL = preReads + postReads) %>% arrange(-FL)

# protein-coding genes
numProteinCodingGenes <- length(unique(class.files$glob_SQ_proteinGenes$associated_gene))
message("Number of protein-coding genes detected:", numProteinCodingGenes)
mostAbundantTranscript <- as.data.frame(class.files$glob_SQ_proteinGenes) %>% 
  # select transcript with the highest number of nreads and supported by nsamples
  group_by(associated_gene) %>%
  arrange(desc(whole_nreads), desc(whole_nsamples)) %>% dplyr::slice(1) %>% ungroup()
mostAbundantTranscript <- as.data.frame(mostAbundantTranscript)
mostAbundantTranscript %>% group_by(structural_category) %>% tally(n = "numTranscripts") %>% mutate(perc = numTranscripts/sum(numTranscripts) * 100)
numProteinCodingGenesFSMTranscript <- length(unique(mostAbundantTranscript[mostAbundantTranscript$structural_category == "FSM","associated_gene"]))
message("Number of protein-coding genes with FSM as most abundant transcript: ", numProteinCodingGenesFSMTranscript)
message("% of protein-coding genes with FSM as most abundant transcript: ", numProteinCodingGenesFSMTranscript/numProteinCodingGenes * 100)

# only select for the base identifier (ignore the decimal point as that's for different versions of the same transcript)
mostAbundantTranscriptFSM <- mostAbundantTranscript %>% filter(structural_category == "FSM") %>% 
  mutate(transcript_firstPart = word(associated_transcript,c(1), sep = fixed(".")))
MANEselectGenes <- MANEselectGenes %>% mutate(Ensembl_nuc_firstPart = word(Ensembl_nuc, c(1), sep = fixed(".")))
nrow(mostAbundantTranscriptFSM %>% filter(associated_transcript %in% MANEselectGenes$Ensembl_nuc))
nrow(mostAbundantTranscriptFSM %>% filter(transcript_firstPart %in% MANEselectGenes$Ensembl_nuc_firstPart))
message("Missing protein-coding genes without a known MALE select:", length(setdiff(mostAbundantTranscriptFSM$associated_gene, MANEselectGenes$symbol)))

# known genes with no FSM transcripts
GenesWithFSM <- unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "FSM","associated_gene"])
GeneswithnoFSM <- setdiff(unique(class.files$glob_SQ_annoGene$associated_gene),GenesWithFSM)
message("Number of genes with no FSM transcripts:", length(GeneswithnoFSM))
message("Proporition to all known genes: ", length(GeneswithnoFSM)/ length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100)

# protein coding genes with more than one transcript
ProteinCodingGenes1 <- class.files$glob_SQ_proteinGenes %>% group_by(associated_gene) %>% tally() %>% filter(n == 1)
ProteinCodingGenesMoreThan1 <- class.files$glob_SQ_proteinGenes %>% group_by(associated_gene) %>% tally() %>% filter(n > 1)
class.files$glob_SQ_proteinGenes_morethan1 <- class.files$glob_SQ_proteinGenes %>% filter(associated_gene %in% ProteinCodingGenesMoreThan1$associated_gene)

# protein coding genes with more than one transcript and no FSM transcripts
ProteinGenesWithFSM <- unique(class.files$glob_SQ_proteinGenes_morethan1[class.files$glob_SQ_proteinGenes_morethan1$structural_category == "FSM","associated_gene"])
ProteinGeneswithnoFSM <- setdiff(unique(class.files$glob_SQ_proteinGenes_morethan1$associated_gene),ProteinGenesWithFSM)

message("Number of protein-coding genes more than one transcript but no FSM transcripts:", length(ProteinGeneswithnoFSM))
message("Proportion to all protein-coding genes: ", length(ProteinGeneswithnoFSM)/ length(unique(class.files$glob_SQ_proteinGenes$associated_gene)) * 100)

class.files$glob_SQ_annoGene %>% filter(structural_category %in% c("NIC","NNC")) %>% arrange(-whole_nreads)
# ONT2.3331.23589 <- ONT2_3223_23974
# ONT2.3331.23546 <- ONT2_3223_23931

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

message("Number of transcripts in prenatal with more than 50 reads: ", length(unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 50 & class.files$glob_SQ_annoGene$postReads == 0,"isoform"])))
message("Number of transcripts in postnatal with more than 50 reads: ", length(unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 50 & class.files$glob_SQ_annoGene$preReads == 0,"isoform"])))

novelTranscripts <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel",]
message("Number of novel transcripts detected in both prenatal and postnatal: ", nrow(novelTranscripts[novelTranscripts$DevStatus == "Both",]))
message("%: ", nrow(novelTranscripts[novelTranscripts$DevStatus == "Both",])/nrow(novelTranscripts) * 100)

genesPrenatal <- unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$DevStatus == "prenatal","associated_gene"])
genesPrePostnatal <- unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$DevStatus == "Both","associated_gene"])
# add RP11-1I2.1 with one transcript
genesPrePostnatal <- c(genesPrePostnatal, "RP11-1I2.1")
genesPostnatal <- unique(class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$DevStatus == "postnatal","associated_gene"])
genesPrenatalOnly <- setdiff(genesPrenatal, c(genesPrePostnatal, genesPostnatal))
genesPostnatalOnly <- setdiff(genesPostnatal,c(genesPrePostnatal, genesPrenatal))

message("Genes uniquely expressed in prenatal: ", length(genesPrenatalOnly))
message("Genes unique expressed in prenatal and protein-coding: ", length(intersect(genesPrenatalOnly, unique(class.files$glob_SQ_proteinGenes$associated_gene))))
message("Proportion: ", length(genesPrenatalOnly)/length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100)

message("Genes uniquely expressed in postnatal: ", length(genesPostnatalOnly))
message("Genes unique expressed in postnatal and protein-coding: ", length(intersect(genesPostnatalOnly, unique(class.files$glob_SQ_proteinGenes$associated_gene))))
message("Proportion: ", length(genesPostnatalOnly)/length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100)

length(genesPrePostnatal)
length(genesPrePostnatal)/length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100

# santiy check of genes capture correctly above
length(genesPrenatalOnly) + length(genesPostnatalOnly) + length(genesPrePostnatal) == length(unique(class.files$glob_SQ_annoGene$associated_gene))
setdiff(unique(class.files$glob_SQ_annoGene$associated_gene), c(genesPrenatalOnly, genesPostnatalOnly, genesPrePostnatal))
genes <- c(genesPostnatalOnly, genesPrenatalOnly, genesPrePostnatal)
setdiff(unique(class.files$glob_SQ_annoGene$associated_gene), genes)
setdiff(genes, unique(class.files$glob_SQ_annoGene$associated_gene))

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

# ERCC
ercc_usage <- plot_usage_persample(ercc.class.file, ercc.demux)

# fusion
message("Number of fusion transcripts: ", nrow(class.files$glob_SQ[class.files$glob_SQ$structural_category == "Fusion",]))
message("% of fusion transcripts: ", nrow(class.files$glob_SQ[class.files$glob_SQ$structural_category == "Fusion",])/nrow(class.files$glob_SQ) * 100)
message("Number of genes involved: ", length(unique(class.files$glob_SQ[class.files$glob_SQ$structural_category == "Fusion","associated_gene"])) * 2)
message("% of genes involved: ", length(unique(class.files$glob_SQ[class.files$glob_SQ$structural_category == "Fusion","associated_gene"])) * 2/ length(unique(class.files$glob_SQ_annoGene$associated_gene)) * 100)
message("number of transcripts from stringent filter fusion: ", nrow(fusion))
message("Number of genes involved: ", length(unique(c(fusion$gene1, fusion$gene2))))

## ------- whole vs targeted dataset ---------

# matchSumTargeted generated from function
whole_vs_targeted_plots(classfiles=class.files$glob_targ_SQ, wholeSamples=wholematchedsamples, targetedSamples=manifest[manifest$ID %in% targetedmatchedsamples,"Sample"], targetGene=selectedTargetGenes)
message("Number of transcripts detected whole only: ", length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Whole","isoform"])))
message("Number of transcripts detected Targeted only: ", length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Targeted","isoform"])))
message("Number of transcripts detected in both whole and targeted only: ", length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Both","isoform"])))

message("Number of novel transcripts detected in Targeted only: ", length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Targeted" & !matchedSumTargeted$structural_category %in% c("FSM","ISM"),"isoform"])))
              
intersect(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Targeted","isoform"]),
          unique(matchedSumTargeted[matchedSumTargeted$dataset == "Both","isoform"]))
intersect(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Whole","isoform"]),
          unique(matchedSumTargeted[matchedSumTargeted$dataset == "Both","isoform"]))

message("Percentage of transcripts validated in whole: ", 
        round(length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Both","isoform"]))/
          length(unique(matchedSumTargeted[matchedSumTargeted$dataset %in% c("Both","Whole"),"isoform"])) * 100,2))
message("Percentage of novel transcripts validated in whole * 100: ",
length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Both" & !matchedSumTargeted$structural_category %in% c("FSM","ISM"),"isoform"]))/
length(unique(matchedSumTargeted[matchedSumTargeted$dataset %in% c("Both","Whole") & !matchedSumTargeted$structural_category %in% c("FSM","ISM"),"isoform"])))
length(unique(matchedSumTargeted[matchedSumTargeted$dataset == "Targeted" & !matchedSumTargeted$structural_category %in% c("FSM","ISM"),"isoform"]))


## ------- long-read sequencing datasets comparisons -------

# number of overlaps between Leung et al. 2021 and whole dataset

comparison_dataset <- function(gfftmap, classfiles, altDataset){
  
        # keep only those that are in the 10reads, 10samples classfiles
        gfftmap <- gfftmap %>% filter(ref_id %in% classfiles$isoform)

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
        
        output <- list(Unique_RB, Common_RB)
        names(output) <- c("Unique_RB","Common_RB")
        return(output)

}

cellReportFurther <- comparison_dataset(gfftmapComparisons$cellReports, class.files$glob_SQ, humanCTX)
directRNAFurther <- comparison_dataset(gfftmapComparisons$directRNA, class.files$glob_SQ, directRNA)
comparison_dataset(gfftmapComparisons$BDRNatureComms, class.files$glob_targ_SQ_20AD, BDRNatureComms)
comparison_dataset(gfftmapComparisons$Herberle, class.files$glob_SQ, HerberleCTX)
comparison_dataset(gfftmapComparisons$Patowary, class.files$glob_SQ, PatowaryCTX)
comparison_dataset(gfftmapComparisons$PatowaryHerbele, PatowaryCTX, HerberleCTX)

# expression difference between CellReports and Bamford dataset for unique transcripts
ExpressionDiffPacBio <- rbind(cellReportFurther$Unique_RB %>% mutate(dataset = "Unique") %>% select(whole_nreads, dataset),
                              cellReportFurther$Common_RB %>% mutate(dataset = "Common") %>% select(whole_nreads, dataset))
ExpressionDiffPacBio %>% group_by(dataset) %>% summarise(median = median(whole_nreads))
res <- wilcox.test(whole_nreads ~ dataset, ExpressionDiffPacBio, exact = FALSE)

res$p.value
format(res$p.value, scientific = TRUE)


# expression difference between cDNA and dRNA for unique transcripts
ExpressionDiffcDNAdRNA <- rbind(directRNAFurther$Unique_RB %>% mutate(dataset = "Unique") %>% select(whole_nreads, dataset),
                                directRNAFurther$Common_RB %>% mutate(dataset = "Common") %>% select(whole_nreads, dataset))
ExpressionDiffcDNAdRNA %>% group_by(dataset) %>% summarise(median = median(whole_nreads))
res <- wilcox.test(whole_nreads ~ dataset, ExpressionDiffcDNAdRNA, exact = FALSE)
res
res$p.value
format(res$p.value, scientific = TRUE)

# novel genes and fusion transcripts with direct RNA validation 
fusion <- fusion %>% mutate(directRNA = ifelse(isoform %in% directRNAFurther$Common_RB$isoform, TRUE,FALSE))
intergenic <- intergenic %>% mutate(directRNA = ifelse(isoform %in% directRNAFurther$Common_RB$isoform, TRUE,FALSE))

write.csv(intergenic, paste0(output_dir, "intergenic.csv"))
write.csv(fusion, paste0(output_dir, "fusion.csv"))

## ------- differential transcript expression -------

# development
nrow(WholeDTE$age) == length(unique(WholeDTE$age$isoform))
message("Number of differentially expressed isoforms: ", nrow(WholeDTE$age))
message("Number of genes with differentially expressed isoforms: ", length(unique(WholeDTE$age$associated_gene)))
message("Number of upregulated transcripts in postnatal vs prenatal: ", nrow(WholeDTE$age[WholeDTE$age$dirAcrossDev == "upregulated",]))
message("top-ranked transcript differentially expressed across development")
head(WholeDTE$age %>% arrange(padj))

# number of upregulated
message("Number of upregulated isoforms: ",nrow(WholeDTE$age %>% filter(dirAcrossDev == "upregulated")))
message("Proportion upregulated isoforms: ",nrow(WholeDTE$age %>% filter(dirAcrossDev == "upregulated"))/nrow(WholeDTE$age))

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


## ---------- Differential transcript usage ----------

head(DIUSig$wholeAge %>% filter(DGE_Dev == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, -as.numeric(totalChange)))
message("Number of genes with significant DIU across development: ", nrow(DIUSig$wholeAge))
message("Number of genes with significant DIU across development and podium Change: ", nrow(DIUSig$wholeAge %>% filter(podiumChange == TRUE)))
message("Number of genes with significant DIU across development and not podium Change: ", nrow(DIUSig$wholeAge %>% filter(podiumChange == FALSE)))
message("Number of genes with significant DIU across development and not DGE: ", nrow(DIUSig$wholeAge %>% filter(DGE_Dev == FALSE)))

# GMP6A
GMP6AStats <- plotIFWholebyGene("GMP6A",DIUnormExp$GMP6A_group,facetTranscriptsFeature=TRUE,sexFeature=FALSE)
WholeDESeqGene$age[WholeDESeqGene$age$associated_gene == "GPM6A",]
t.test(perc ~ group, GPM6AStats[GPM6AStats$isoform == "ONT4.13313.18545",])
t.test(perc ~ group, GPM6AStats[GPM6AStats$isoform == "ONT4.13313.3469",])
t.test(perc ~ group, GPM6AStats[GPM6AStats$isoform == "ONT4.13313.3457",])
t.test(perc ~ group, GPM6AStats[GPM6AStats$isoform == "ONT4.13313.3461",])
t.test(perc ~ group, GPM6AStats[GPM6AStats$isoform == "ONT4.13313.18565",])

## ---------- Alternative splicing events ----------

# calculate AS events for pre-natal and post-natal detected transcripts 
splicingExp <- finalTranscriptClassificationTranscript %>% mutate(A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3),
                                                          ES = as.numeric(ES), ESpre = ES * normPre/sum(ES), ESpost = ES * normPost/sum(ES),
                                                          NE = as.numeric(NE_Int), NEpre = NE_Int * normPre/sum(NE_Int), NEpost = NE * normPost/sum(NE_Int),
                                                          IR = as.numeric(IR), IRpre = IR * normPre/sum(IR), IRpost = IR * normPost/sum(IR),
                                                          A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3))

# subset events
FICLEES <- splicingExp %>% select(isoform, ESpre, ESpost) %>% reshape2::melt(variable.name = "development", value.name = "normES") 
FICLEIR <- splicingExp %>% select(isoform, IRpre, IRpost) %>% reshape2::melt(variable.name = "development", value.name = "normIR") 
FICLEA5A3 <- splicingExp %>% select(isoform, A5A3pre, A5A3post) %>% reshape2::melt(variable.name = "development", value.name = "normA5A3") 
FICLENE <- splicingExp %>% select(isoform, NEpre, NEpost) %>% reshape2::melt(variable.name = "development", value.name = "normNE") 

ggplot(FICLEES, aes(x = development, y = log10(normES))) + geom_boxplot()
ggplot(FICLEIR, aes(x = development, y = log10(normIR))) + geom_boxplot()
ggplot(FICLEA5A3, aes(x = development, y = log10(normA5A3))) + geom_boxplot()

t.test(normES ~ development, FICLEES)
t.test(normA5A3 ~ development, FICLEA5A3)
t.test(normIR ~ development, FICLEIR)
t.test(normNE ~ development, FICLENE)

## ---- exon skipping events

# identify transcripts with exon skipping events (ES > 0)
ES <- finalTranscriptClassificationTranscript[finalTranscriptClassificationTranscript$ES > 0,]

# identify transcripts in whole transcriptome dataset with exon skipping events
class.files$glob_SQ_proteinGenes_ES <- class.files$glob_SQ_proteinGenes[class.files$glob_SQ_proteinGenes$isoform %in% ES$isoform,] 
ES_glob_SQ <- ES %>% filter(isoform %in% class.files$glob_SQ_proteinGenes_ES$isoform)

message("Number of transcripts from protein-coding genes with exon skipping events", nrow(class.files$glob_SQ_proteinGenes_ES))
message("Number of those transcripts by development: ")
table(class.files$glob_SQ_proteinGenes_ES$DevStatus)

message("Total number of exon skipping events", sum(ES_glob_SQ$ES))
message("Number of exon skipping events in prenatal: ", sum(ES_glob_SQ[ES_glob_SQ$DevStatus == "prenatal","ES"]))
message("Number of exon skipping events in postnatal: ", sum(ES_glob_SQ[ES_glob_SQ$DevStatus == "postnatal","ES"]))
message("Number of exon skipping events in postnatal and postnatal: ", sum(ES_glob_SQ[ES_glob_SQ$DevStatus == "Both","ES"]))

message("Median number of exon skipping events per transcript: ", median(ES_glob_SQ$ES))
message("Max number of exon skipping events in any given transcript: ", max(ES_glob_SQ$ES))
message("Top-ranked transcripts with the most exon skipping events: ")
ES_glob_SQ %>% arrange(-ES) %>% head(.)

## ---- novel exon inclusion events

# identify transcripts with novel exons 
# retain only NIC and NNC transcripts 
NE <- finalTranscriptClassificationTranscript[finalTranscriptClassificationTranscript$NE_Int > 0,]
class.files$glob_SQ_proteinGenes_NE <- class.files$glob_SQ_proteinGenes[class.files$glob_SQ_proteinGenes$isoform %in% NE$isoform,] %>% filter(structural_category %in% c("NIC","NNC"))

message("Number of transcripts from protein-coding genes with exon skipping events", nrow(class.files$glob_SQ_proteinGenes_NE))
message("% of all transcripts annotated to protein-coding genes", nrow(class.files$glob_SQ_proteinGenes_NE)/nrow(class.files$glob_SQ_proteinGenes))

message("Number of those transcripts by development: ")
table(class.files$glob_SQ_proteinGenes_NE$DevStatus)

message("Median number of reads (expression) of transcripts with novel exons: ", median(class.files$glob_SQ_proteinGenes_NE$whole_nreads))
message("Top-ranked most abundantly expressed transcripts with novel exons")
class.files$glob_SQ_proteinGenes_NE %>% select(isoform, associated_gene, whole_nreads, whole_nsamples) %>% arrange(-whole_nreads) %>% head(.)


## ---------- Proteogenomics ----------

# high-level coding potential (cpat)
cpat_glob_SQ <- merge(cpat, class.files$glob_SQ[,c("isoform","structural_category")], by.x = "seq_ID", by.y = "isoform")
cpat_glob_SQ <- cpat_glob_SQ  %>% filter(seq_ID %in% class.files$glob_SQ$isoform)
cpat_glob_SQ_high_coding <- cpat_glob_SQ %>% filter(Coding_prob >= 0.364)
message("Number of transcripts with predicted CPAT score: ", length(unique(cpat_glob_SQ$seq_ID)))
message("% of transcripts with predicted CPAT score: ", length(unique(cpat_glob_SQ$seq_ID))/nrow(class.files$glob_SQ))
message("Number of transcripts with high prediction score: ", length(unique(cpat_glob_SQ_high_coding$seq_ID)))
message("% of transcripts with high prediction score, NIC: ", 
        nrow(cpat_glob_SQ_high_coding[cpat_glob_SQ_high_coding$structural_category == "NIC",])/length(unique(cpat_glob_SQ$seq_ID)) * 100)
message("% of transcripts with high prediction score, NNC: ", 
        nrow(cpat_glob_SQ_high_coding[cpat_glob_SQ_high_coding$structural_category == "NNC",])/length(unique(cpat_glob_SQ$seq_ID)) * 100)

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

novelPeptides <- novelPeptides %>% mutate(length = nchar(as.character(seq)))
novelPeptidesunique <- novelPeptides %>% group_by(acc) %>% top_n(1, length)
colnames(novelPeptides) <- c("gene","isoform","protein_seq")
novelPeptides <- merge(novelPeptides,class.files$glob_SQ, by = "isoform") %>% select(isoform, associated_gene, structural_category, protein_seq)
novelPeptidesTranscripts <- novelPeptides[!novelPeptides$structural_category %in% c("FSM","ISM"),]
message("Total number of novel trancripts with mass-spec validation:", length(unique(novelPeptidesTranscripts$isoform)))
write.csv(novelPeptidesTranscripts, paste0(output_dir,"novelTranscriptsPeptides.csv"), quote = F, row.names = F)


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
class.files$glob_SQ[class.files$glob_SQ$associated_transcript == "novel",] %>% arrange(-whole_nreads) %>% select(isoform, structural_category, exons, associated_gene, associated_transcript, nreads, nsamples) %>% head(.)
class.files$bambu[class.files$bambu$associated_gene == "RPS27A",]

# isoform fraction of novel transcript
tabulateIF(class.files$glob_SQ %>% filter(associated_gene %in% "RPS27A"), "Whole") %>% 
  mutate(label = ifelse(perc > 10, isoform, NA), perc = perc / 100)  %>% filter(isoform == "ONT2.3331.23589")
tabulateIF(class.files$glob_SQ %>% filter(associated_gene %in% "RPS27A"), "Whole") %>% 
  mutate(label = ifelse(perc > 10, isoform, NA), perc = perc / 100)  %>% filter(isoform == "ONT2.3331.23546")


novelPeptidesTranscripts <- read.table(paste0(dirnames$output,"novelTranscriptsPeptides.txt"))
class.files$glob_targ_SQ_counts[class.files$glob_targ_SQ_counts$isoform %in% novelPeptidesTranscripts$isoform,] %>% arrange(-nreads) 
class.files$bambu[class.files$bambu$associated_gene == "HMOX2",]
  
