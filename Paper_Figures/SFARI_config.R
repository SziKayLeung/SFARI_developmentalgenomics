#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("vroom"))

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN,"merge_characterise_dataset/run_ggtranscript.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)


## ------------ directory names --------------- 

root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
root_rb_dir <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/"
recovered_dir <- "/lustre/recovered/Research_Project-MRC148213/sl693/RBFetal/"
dirnames <- list(
  
  wholetarg_SQ = paste0(root_rb_dir,"6_sqanti3/"),
  
  output = paste0(root_sfari,"/0_output/"),
  utils = paste0(root_sfari,"/0_utils/"),
  protein = paste0(root_sfari, "/8_longReadProteogenomics/longReadProteogenomics"),
  
  DGE = paste0(root_sfari, "10_deseq//1_DGE/"), 
  DTE = paste0(root_sfari, "10_deseq//2_DTE/"),
  DIU = paste0(root_sfari, "10_deseq//3_DIU/"),
  
  ficle = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/TargetGenes",
  
  # Leung et al. 2021 PacBio HumanCTX dataset
  humanPacBio = "/lustre/projects/Research_Project-MRC148213/lsl693/PacBioPaper/SQANTI2/HumanCTX"
)

TargetGene = read.table(paste0(root_sfari, "0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]
# list of SZ and ASD genes in targeted panel 
#TargetGeneSZASD = read.csv(paste0(root_sfari,"0_metadata/TargetGeneByDisease.csv"))


## ------------- Phenotype files -------------------
phenotype <- fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
phenotype <- phenotype %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), col = paste0(sample,"_",group)) 
femaleWhole <- phenotype[phenotype$sex == "F" & grepl("Whole",phenotype$sample),][["sample"]]
maleWhole <- phenotype[phenotype$sex == "M" & grepl("Whole",phenotype$sample),][["sample"]]
postWhole <- phenotype[phenotype$group == "Postnatal" & grepl("Whole",phenotype$sample),][["sample"]]
preWhole <- phenotype[phenotype$group == "Prenatal" & grepl("Whole",phenotype$sample),][["sample"]]


## -------------- Final classification files ------------- 

# class.files
# list: glob_targ_SQ, glob_SQ, targ_SQ
load(file = paste0(dirnames$wholetarg_SQ,"all_filtered_classification_2reads2samples_noMonoIntergenicAll.RData"))
class.files$protein_filtered = fread(paste0(dirnames$protein,"/7_classified_protein/Whole.sqanti_protein_classification.tsv"),data.table = F)

# subset by number of FL reads
femaleReads <- class.files$glob_SQ %>% select(all_of(femaleWhole)) %>% apply(., 1, sum)
maleReads <- class.files$glob_SQ %>% select(all_of(maleWhole)) %>% apply(., 1, sum)
postReads <- class.files$glob_SQ %>% select(all_of(postWhole)) %>% apply(., 1, sum)
preReads <- class.files$glob_SQ %>% select(all_of(preWhole)) %>% apply(., 1, sum)
class.files$glob_SQ <- class.files$glob_SQ %>% mutate(FReads = femaleReads, MReads = maleReads, preReads = preReads, postReads = postReads)
class.files$glob_SQ$DevStatus <- apply(class.files$glob_SQ, 1, function(x) identify_dataset_by_counts(x[["postReads"]], x[["preReads"]], "postnatal","prenatal"))

# annotated genes 
class.files$glob_SQ_annoGene <- class.files$glob_SQ %>% filter(!grepl("novelGene", associated_gene))
class.files$glob_SQ_annoGene <- class.files$glob_SQ_annoGene %>% mutate(novelTranscript = ifelse(associated_transcript == "novel","Novel","Known"))

wholesamples <- colnames(class.files$glob_targ_SQ)[grepl("Whole", colnames(class.files$glob_targ_SQ))]
targetedsamples <- colnames(class.files$glob_targ_SQ)[grepl("Targeted", colnames(class.files$glob_targ_SQ))]
matchedsamples <- intersect(gsub("^.*?Whole","",wholesamples),gsub("^.*?Targeted","",targetedsamples))
targetedmatchedsamples <- paste0("Targeted",matchedsamples)
wholematchedsamples <- paste0("Whole",matchedsamples)


## -------------- DESeq2 ----------------

# WholeDESeq
WholeDTE <- list(
  sex = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_sex_resSig.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_development_resSig.csv"),delim = ",",show_col_types = FALSE)
)
WholeDTE <- lapply(WholeDTE, function(x) x %>% mutate(dirAcrossDev = ifelse(log2FoldChange < 0 , "upregulated", "downregulated")))

normWhole <- fread(paste0(dirnames$DTE,"DESeq2_whole_development_normAll.csv"))
WholePreNorm <- normWhole %>% filter(sample %in% preWhole) %>% group_by(isoform) %>% tally(normalised_counts)
WholePostNorm <- normWhole %>% filter(sample %in% postWhole) %>% group_by(isoform) %>% tally(normalised_counts)
normWholeIsoform <- merge(WholePreNorm %>% `colnames<-`(c("isoform", "normPre")) , WholePostNorm %>% `colnames<-`(c("isoform", "normPost")), by = "isoform")

## -------------- differenetial gene expression -------------
WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(dirnames$DGE,"DESeq2_whole_gene_sex_resSig.csv"))),
  age = as.data.frame(fread(paste0(dirnames$DGE,"DESeq2_whole_gene_development_resSig.csv")))
)


## -------------- differential isoform usage ----------------

load(file = paste0(dirnames$output,"DIUSig.RData"))
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
                                                    DGE_Sex = ifelse(Gene %in% WholeDESeqGeneSig$sex$associated_gene,TRUE,FALSE))


## -------------- normalized counts ----------------

#Exp
load(file = paste0(dirnames$DTE,"DESeq2_whole_normSig.RData"))
load(file = paste0(dirnames$DTE,"DESeq2_whole_normAll.RData"))
load(file = paste0(dirnames$DGE,"DESeq2_whole_normAll.RData"))

# raw counts 
rawCounts <- read.table(paste0(root_sfari, "/4_transcriptClean/cluster_report_counts.txt"), col.names = c("counts","file")) 
rawCounts <- rawCounts %>% mutate(sample = word(file,c(1),sep=fixed(".")))
# corrected phenotype sample to match
rawCounts$sample[rawCounts$sample == "WholeN51"] <- "WholeN55"
rawCounts <- merge(rawCounts, phenotype, by = "sample", all = T)

## -------------- cpat ----------------

Cpat <- list(
  whole = data.table::fread(paste0(dirnames$protein,"/5_calledOrfs/Whole.ORF_prob.best.tsv"))#,
  #wholeTargeted = data.table::fread(paste0(recovered_dir,"/2_cpat_tc20bp/WholeTargeted.ORF_prob.best.tsv"))
)
#save(Cpat, file = paste0(root_dir,"RBFetal/2_cpat_tc20bp/Whole.ORF_prob.best.RData"))
# keep  only the list of isoforms in the final dataset
Cpat$whole <- Cpat$whole %>% filter(seq_ID %in% class.files$glob_SQ$isoform) %>% mutate(coding_status = ifelse(Coding_prob >= 0.364, "Coding","Non_Coding"))
Cpat$whole_noORF <- read.table(paste0(dirnames$protein,"/5_calledOrfs/Whole.no_ORF.txt")) %>% mutate(coding_status = "No_ORF") %>% `colnames<-`(c("seq_ID", "coding_status"))

## -------------- gtf ----------------

gtf <- list(
  glob_targ = rtracklayer::import(paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf"))
  #ref = rtracklayer::import(paste0(dirnames$output,"gencode.v40.ggTransRefGenes.gtf"))
  #ref = rtracklayer::import(paste0(dirnames$output,"refExons.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))
gtf$ref <- data.table::fread(paste0(dirnames$utils,"refExons.gtf")) %>% dplyr::rename("gene_id" = "gene_name") %>% mutate(type = "exon")
gtf$merged <- rbind(gtf$glob_targ[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                         gtf$ref[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

GI <- c("GRIN2A","GRIA3","SEPTIN4","RTN4","MBP","RPS4Y1","XIST","ADD3","CNTNAP2","ANKRD12","VXN","PKM")
RefIsoforms <- lapply(GI, function(x) unique(gtf$ref[gtf$ref$gene_id == x & !is.na(gtf$ref$transcript_id), "transcript_id"]))
names(RefIsoforms ) <- GI


## -------------- disease list ----------------

diseasegenelists <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/0_utilities/disease_list/"
SFARI <- read.csv(paste0(diseasegenelists,"SFARI-Gene_genes_07-17-2023release_09-26-2023export.csv"),header=T)
SFARI_CLASS_1_2_S <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2 |SFARI$syndromic == 1)
SFARI_CLASS_1_2 <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2)
disease_list <- list(
  SCHEMA = read.table(paste0(diseasegenelists,"SCHEMA_Oct2023.csv"), sep=",", header=T, stringsAsFactors = F),
  DDG2P = read.table(paste0(diseasegenelists,"DDG2P_26_9_2023.csv"), sep=",", header=T, stringsAsFactors = F) %>% filter(confidence.category != "limited")
)


# X inactivation list 
XExcapeList <- read.csv(paste0(dirnames$utils,"XInactivation.csv"))


## -------------- long read proteoogenomics ----------------

load(paste0(dirnames$utils,"/proteinInputWhole.RData"))

# filter protein input to only the class files (pb_accs as this was the original pb_accs before being collapsed)
proteinInput$t2p.collapse <- proteinInput$t2p.collapse %>% filter(pb_accs %in% class.files$glob_SQ$isoform)

proteinInput$filtered <- read.table(paste0(dirnames$protein,"/7_classified_protein/Whole.classification_filtered.tsv"), sep = "\t", header = T)

ProteinInput = list(
  cpat = read.table(paste0(dirnames$mprotein,"5_calledOrfs/all_iso_ont.ORF_prob.best.tsv"), sep ="\t", header = T),
  cpat_best = read.table(paste0(dirnames$mprotein,"5_calledOrfs/all_iso_ont_best_orf.tsv"), sep ="\t", header = T),
  mapped = read.table(paste0(dirnames$mprotein,"5_calledOrfs/all_orfs_mapped.tsv"), sep ="\t", header = T),
  noORF = read.table(paste0(dirnames$mprotein,"5_calledOrfs/all_iso_ont.no_ORF.txt"), sep ="\t", header = F),
  t2p.collapse = read.table(paste0(dirnames$mprotein,"6_refined_database/all_iso_ont_orf_refined.tsv"), sep = "\t", header = T),
  t2p.collapse.refined = read.table(paste0(dirnames$mprotein,"6_refined_database/all_iso_ont_orf_refined_collapsed.tsv"))
)



## -------------- prenatal vs postnatal ----------------

class.files$glob_SQ_annoGene_prenatal <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,]
class.files$glob_SQ_annoGene_postnatal <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,]


## -------------- comparison with Leung et al.(2021) dataset ----------------

# gffcompare output 
gfftmap <- data.table::fread(paste0(root_sfari,"14_OverlapPacBio/sfariPacBio.HumanCTX.collapsed_classification.filtered_lite.gtf.tmap"), data.table = FALSE)
humanCTX <- read.table(paste0(dirnames$humanPacBio, "/HumanCTX.collapsed_classification.filtered_lite_classification.txt"), header = TRUE)
humanCTX$totalFL <- humanCTX %>% dplyr::select(contains("FL.")) %>% apply(.,1,sum)


## -------------- comparison with Leung et al.(2021) dataset ----------------

protein_coding_genes = read.table("/lustre/home/vc362/protein-coding-genes.txt")
FICLE_class <- fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/TargetGenes/all_final_transcript_classifications.csv", data.table = F)
FICLE_class <- FICLE_class[!duplicated(FICLE_class), ]
FICLE_class <- FICLE_class %>% select(-isoform)
FICLE_class_final <- FICLE_class %>% filter(isoform %in% class.files$glob_SQ_annoGene$isoform)
protein_coding_genes_isoforms <- class.files$glob_SQ[class.files$glob_SQ$associated_gene %in% protein_coding_genes$V1,"isoform"]
setdiff(rownames(FICLE_class_final$isoform), protein_coding_genes_isoforms)
FICLENotProcessed <- setdiff(protein_coding_genes_isoforms, rownames(FICLE_class_final$isoform))
FICLEProcessedGenes <- unique(class.files$glob_SQ[class.files$glob_SQ$isoform %in% FICLE_class_final$isoform,"associated_gene"])
FICLENotProcessedGenes <- unique(class.files$glob_SQ[class.files$glob_SQ$isoform %in% FICLENotProcessed,"associated_gene"])
unique(setdiff(FICLENotProcessedGenes,FICLEProcessedGenes))

FICLE_class_final_exp <- merge(FICLE_class_final[,c("isoform","A5A3","ES","IR","NE_All")], class.files$glob_SQ[,c("isoform","preReads","postReads")], by = "isoform")
FICLE_class_final_exp <- merge(FICLE_class_final_exp, normWholeIsoform, by = "isoform", all.x = T)
