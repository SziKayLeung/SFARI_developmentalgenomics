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
LOGEN_ROOT <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGEN,"compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN,"merge_characterise_dataset/run_ggtranscript.R"))
sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"longread_QC"), pattern="*.R", full = T), source,.GlobalEnv)
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)


## ------------ directory names --------------- 

root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
root_rna_dir <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/"
dirnames <- list(
  
  # general
  output = paste0(root_sfari,"/0_output/"),
  utils = paste0(root_sfari,"/0_utils/"),

  wholetarg_SQ = paste0(root_sfari,"C_Whole_Targeted/9_sqanti_final/"),
  protein = paste0(root_sfari, "C_Whole_Targeted/10_longReadProteogenomics/"),
  massspec= paste0(root_sfari, "C_Whole_Targeted/13_massSpec/"),
  
  # whole dataset differential expression analysis
  DGE = paste0(root_sfari, "A_Whole/10_deseq/1_DGE/"), 
  DTE = paste0(root_sfari, "A_Whole/10_deseq/2_DTE/"),
  DIU = paste0(root_sfari, "A_Whole/10_deseq/3_DIU/"),

  # ficle
  ficle = paste0(root_sfari, "C_Whole_Targeted/15_ficle/TargetGenes/"),

  # comparison to other long read sequencing datasets
  humanPacBio = paste0(root_sfari, "C_Whole_Targeted/17_longReadDatasetComparisons/Leung2021/HumanCTX"),
  directRNA = paste0(root_rna_dir, "dRNA/Rosie/9_sqanti_final/"),

  overlapDatasets = paste0(root_sfari, "C_Whole_Targeted/14_OverlapDatasets/")
)

TargetGene = read.table(paste0(root_sfari, "0_metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]
ProteinCodingGenes = read.table(paste0(dirnames$utils, "WholeProteinCodingGenes.txt"))[["V1"]]
GWAS = c("ACTR1B", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CSMD1", "CTD-2008L17.2", "CUL9", "DCC", "DLGAP2", "DPYD", "EMX1", "ENOX1", "EPN2", "EYS", "FURIN", "GABBR2", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "IL1RAPL1", "IMMP2L", "IRF3", "KIAA1549", "KLF6", "LINC00320", "LINC01088", "LRRC4B", "MAD1L1", "MAN2A1", "MAPT", "MSI2", "NAB2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCGF3", "PCNXL3", "PDE4B", "PJA1", "PLCH2", "PTPRD", "R3HDM2", "RP11-399D6.2", "RP11-507B12.2", "SGCD", "SLC39A8", "SLC4A10", "SNAP91", "SP4", "THAP8", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WSCD2", "ZNF804A", "ZNF823", "ZNF835")


## ------------- Phenotype files -------------------

phenotype <- fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
phenotype <- phenotype %>% mutate(type = ifelse(grepl("Targeted",sample),"Targeted","Whole"),
                     sampleID = gsub("^Targeted", "", sample)) %>% mutate(sampleID = gsub("^Whole","", sampleID))
phenotype <- phenotype %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), 
                                  col = paste0(sample,"_",group),
                                  sex = factor(sex, levels = c("M","F"))) 
femaleWhole <- phenotype[phenotype$sex == "F" & grepl("Whole",phenotype$sample),][["sample"]]
maleWhole <- phenotype[phenotype$sex == "M" & grepl("Whole",phenotype$sample),][["sample"]]
postWhole <- phenotype[phenotype$group == "Postnatal" & grepl("Whole",phenotype$sample),][["sample"]]
preWhole <- phenotype[phenotype$group == "Prenatal" & grepl("Whole",phenotype$sample),][["sample"]]
matchedsamples <- intersect(phenotype[phenotype$type == "Whole","sampleID"],phenotype[phenotype$type == "Targeted","sampleID"])
wholematchedsamples <- phenotype[phenotype$sampleID %in% matchedsamples & phenotype$type == "Whole","sample"]
targetedmatchedsamples <- phenotype[phenotype$sampleID %in% matchedsamples & phenotype$type == "Targeted","sample"]

# manifest 
manifest <- fread(paste0(root_sfari, "0_metadata/WholeTargetedphenotype_manifest.csv"),data.table=F, stringsAsFactors=F)

## -------------- Final classification files ------------- 

load(file = paste0(dirnames$wholetarg_SQ,"sqantifiltered_monoexonicfiltered_2reads2samples.RData"))
annoGenesStats <- list(
  novelTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel",],
  annoTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript != "novel",],
  NIC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NIC",],
  NNC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NNC",]
)

## Protein level 
class.files$protein_filtered = fread(paste0(dirnames$protein,"/7_classified_protein/Whole.sqanti_protein_classification.tsv"),data.table = F)

## -------------- Bambu ---------------- 

# bambu collapsed from pbmm2 aligned files (whole + targeted)
#bambu <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/16_bambu/sqanti/WholeTargeted_RulesFilter_result_classification.txt"
#class.files$bambu <- SQANTI_class_preparation(bambu,"nstandard")
#class.files$bambu_annoGene <- class.files$bambu %>% filter(!grepl("novelGene", associated_gene))
#class.files$bambu_annoGene <- class.files$bambu_annoGene  %>% mutate(novelTranscript = ifelse(associated_transcript == "novel","Novel","Known"))

#bambuAnnoGeneStats <- list(
#  novelTrans = class.files$bambu_annoGene[class.files$bambu_annoGene$associated_transcript == "novel",],
#  annoTrans = class.files$bambu_annoGene[class.files$bambu_annoGene$associated_transcript != "novel",],
#  NIC = class.files$bambu_annoGene[class.files$bambu_annoGene$structural_category == "NIC",],
#  NNC = class.files$bambu_annoGene[class.files$bambu_annoGene$structural_category == "NNC",]
#)


## -------------- DESeq2 ----------------

# WholeDESeq
WholeDTE <- list(
  sex = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_sex_resSig.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_development_resSig.csv"),delim = ",",show_col_types = FALSE)
)
WholeDTE <- lapply(WholeDTE, function(x) x %>% mutate(dirAcrossDev = ifelse(log2FoldChange < 0 , "upregulated", "downregulated")))
WholeDTE$sex <- merge(WholeDTE$sex, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)
WholeDTE$age <- merge(WholeDTE$age, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)

WholeDTEAll <- list(
  sex = vroom(paste0(dirnames$DTE,"DESeq2_whole_transcript_sex_resAll.csv"),delim = ",",show_col_types = FALSE)
)
WholeDTEAll$sex <- merge(WholeDTEAll$sex, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)

normWhole <- fread(paste0(dirnames$DTE,"DESeq2_whole_development_normAll.csv"))
WholePreNorm <- normWhole %>% filter(sample %in% preWhole) %>% group_by(isoform) %>% tally(normalised_counts)
WholePostNorm <- normWhole %>% filter(sample %in% postWhole) %>% group_by(isoform) %>% tally(normalised_counts)
normWholeIsoform <- merge(WholePreNorm %>% `colnames<-`(c("isoform", "normPre")) , WholePostNorm %>% `colnames<-`(c("isoform", "normPost")), by = "isoform")

## -------------- differenetial gene expression -------------
WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(dirnames$DGE,"DESeq2_whole_gene_sex_resSig.csv"))),
  age = as.data.frame(fread(paste0(dirnames$DGE,"DESeq2_whole_gene_development_resSig.csv")))
)

#notDGE <- setdiff(class.files$glob_targ_SQ$associated_gene, WholeDESeqGeneSig$age)
#notDGE[!grepl("novelGene",notDGE)]
#write.table(notDGE[!grepl("novelGene",notDGE)],"NotDGEList.csv",row.names=F,col.names = F, sep = ",",quote=F)


## -------------- differential isoform usage ----------------

load(file = paste0(dirnames$output,"DIUSig.RData"))
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
                                                    DGE_Sex = ifelse(Gene %in% WholeDESeqGeneSig$sex$associated_gene,TRUE,FALSE))

DIUSig$wholeAllSex <- DIUSig$wholeAllSex %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
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

GI <- c("GRIN2A","GRIA3","SEPTIN4","RTN4","MBP","RPS4Y1","XIST","ADD3","CNTNAP2","ANKRD12","VXN","PKM","MORF4L2","GNAS","RSP27A")
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
gfftmapComparisons <- list(
  cellReports = data.table::fread(paste0(dirnames$overlapDatasets,"cellReports2021/sfari_PacBio.HumanCTX.collapsed_classification.filtered_lite.gtf.tmap"), data.table = FALSE),
  directRNA = data.table::fread(paste0(dirnames$overlapDatasets,"directRNA/sfari_dRNA.sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf.tmap"), data.table = FALSE),
  BDRNatureComms = data.table::fread(paste0(dirnames$overlapDatasets,"BDRNatureComms2024/sfari_BDR.ontBDR_collapsed.filtered_counts_filtered.gtf.tmap"), data.table = FALSE)
)
humanCTX <- read.table(paste0(dirnames$humanPacBio, "/HumanCTX.collapsed_classification.filtered_lite_classification.txt"), header = TRUE)
humanCTX$totalFL <- humanCTX %>% dplyr::select(contains("FL.")) %>% apply(.,1,sum)
directRNA <- read.table(paste0(dirnames$directRNA, "sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt"), header = TRUE, sep = "\t", as.is = T)

## -------------- FICLE ----------------

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


## ----- selected disease genes ------
monoAllelicDDP = c('ADNP','ANK2','ANKRD11','ARID1B','ASH1L','AUTS2','BCL11A','BCL11B','CACNA1C','CACNA1G','CACNA1H','CDH1','CDH2','CDK13','CHD2','CHD8','CLCN3','CNOT1','DLG4','DMD','DYRK1A','EIF2S3','EP300','FMR1','FOXP1','FOXP2','GABBR2','GATA3','GATA4','GATA6','GLMN','GRIA3','GRIK2','GRIN1','GRIN2A','GRIN2B','H1-4','HNF1B','HNF4A','IL1RAPL1','ITPR1','KAT6B','KDM6A','KDM6B','KIRREL3','KMT2D','LMNA','MAGI2','MECP2','MED13L','MEIS2','MNX1','NF1','NFIA','NFIB','NLGN3','NLGN4X','NRXN1','OPA1','PAX6','PBX1','POGZ','POLA1','PTEN','QRICH1','RANBP2','RBFOX1','RPL10','SCN1A','SCN2A','SCN8A','SET','SETD1A','SHANK3','SHH','SLC9A9','SMARCE1','SON','SOX9','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','UBE3A','USP7','WFS1','ZBTB20','ZEB2','ZMYM2')

biallelicDDP = c('ASPM','CACNA1G','CHL1','CISD2','CLCN3','CNTNAP2','CTNNA2','DCC','EIF2AK3','EOMES','GLIS3','GPC6','GRIK2','GRIN1','GRIN2A','GRM1','HADH','HERC1','HPSE2','ITPR1','KDM6B','KIAA1109','LMNA','LRBA','MCPH1','NRXN1','NTNG1','ONECUT1','PDIA6','PDX1','PEX16','PPP1R15B','PTF1A','QARS1','RELN','RFX6','SLC39A8','SLF2','TRMT10A','TSPEAR','UBE3B','VPS13B','WFS1','WWOX','ZBTB16','ZFP57')

GWAS = c("ACTR1B", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CSMD1", "CTD-2008L17.2", "CUL9", "DCC", "DLGAP2", "DPYD", "EMX1", "ENOX1", "EPN2", "EYS", "FURIN", "GABBR2", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "IL1RAPL1", "IMMP2L", "IRF3", "KIAA1549", "KLF6", "LINC00320", "LINC01088", "LRRC4B", "MAD1L1", "MAN2A1", "MAPT", "MSI2", "NAB2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCGF3", "PCNXL3", "PDE4B", "PJA1", "PLCH2", "PTPRD", "R3HDM2", "RP11-399D6.2", "RP11-507B12.2", "SGCD", "SLC39A8", "SLC4A10", "SNAP91", "SP4", "THAP8", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WSCD2", "ZNF804A", "ZNF823", "ZNF835")

SFARI = c('ADNP','AGAP1','ANK2','ANKRD11','ARID1B','ASH1L','ASPM','AUTS2','BCL11A','CACNA1C','CACNA1G','CACNA1H','CADPS','CADPS2','CD38','CDH13','CDH2','CDK13','CELF6','CHD2','CHD8','CLTCL1','CNOT1','CNTN4','CNTN6','CNTNAP2','CSMD1','CTNNA2','CYFIP1','DAGLA','DCC','DISC1','DLG4','DLGAP2','DMD','DPYD','DRD2','DRD3','DYRK1A','ELP4','EP300','FBXO40','FMR1','FOXP1','FOXP2','GABBR2','GPC6','GRIA3','GRIK2','GRIK3','GRIN1','GRIN2A','GRIN2B','H1-4','HERC1','IL1RAPL1','IMMP2L','ITPR1','KAT6B','KATNAL2','KCTD13','KDM6A','KDM6B','KIRREL3','LRBA','MAPT','MCM4','MCPH1','MECP2','MED13L','MEIS2','MET','NEGR1','NF1','NFIA','NFIB','NKX2-2','NLGN1','NLGN2','NLGN3','NLGN4X','NR3C2','NRXN1','NTNG1','NXPH1','OXTR','PARD3B','PAX6','PBX1','PCDH10','PCDH9','PHB','PJA1','POGZ','PRKN','PTEN','PTPRT','QRICH1','RBFOX1','RELN','RPL10','RUNX1T1','SCN1A','SCN2A','SCN8A','SET','SETD1A','SEZ6L2','SHANK3','SLC4A10','SLC9A9','SON','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','TSHZ3','UBE3A','USP7','VPS13B','WWOX','ZBTB16','ZBTB20','ZMYM2','ZNF18','ZNF804A')

schemaGenes <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/schema_genes.txt", col.names = c("Gene"))

selectedTargetGenes <- unique(c(as.character(monoAllelicDDP), as.character(SFARI), as.character(schemaGenes$Gene), as.character(biallelicDDP), GWAS))
length(selectedTargetGenes)
