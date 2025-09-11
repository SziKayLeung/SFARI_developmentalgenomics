#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: config file containing variables and paths
##
## ---------------------------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("vroom"))
suppressMessages(library("tidyr"))


## ------------ directory names --------------- 

LOGEN <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/LOGen"
root_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/PaperZenodo/"
output_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/Output/"
source(paste0(LOGEN, "miscellaneous/convert_ensembl_to_symbol.R")) # convert ensembl to gene symbol

TargetGene <- read.table(paste0(root_dir, "metadata/Complete_TargetGenes_TargetedSequencing.txt"))[["V1"]]
TargetGeneS2 <- read.table(paste0(root_dir, "utils/targetGenesSupplementaryTable2.txt"))
ProteinCodingGenes <- read.table(paste0(root_dir, "/utils/protein-coding-genes.txt"))[["V1"]]
message("Known protein coding genes: ", length(ProteinCodingGenes))

MANEselectGenes <- fread(paste0(root_dir, "utils/MANE.GRCh38.v1.4.summary.txt"), data.table = F) %>% 
  filter(MANE_status == "MANE Select")
length(unique(MANEselectGenes$symbol)) == nrow(MANEselectGenes)

# check difference in TargetGeneLists
setdiff(TargetGeneS2$V1,TargetGene)
length(unique(TargetGeneS2$V1))

## ------------- Phenotype files -------------------

phenotype <- fread(paste0(root_dir, "metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
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
manifest <- fread(paste0(root_dir, "metadata/WholeTargetedphenotype_manifest.csv"),data.table=F, stringsAsFactors=F)

## -------------- Final classification files ------------- 

# load glob_targ_SQ, glob_SQ, glob_SQ_annoGene, targ_SQ
load(file = paste0(root_dir,"sqanti/sqantifiltered_monoexonicfiltered_2reads2samples_v2.RData"))
class.files$glob_SQ_annoGene <- class.files$glob_SQ %>% filter(!grepl("novelGene", associated_gene)) %>% mutate(novelTranscript = ifelse(associated_transcript == "novel","Novel","Known"))
class.files$glob_targ_SQ_default_annoGene <- class.files$glob_targ_SQ_default %>% filter(!grepl("novelGene", associated_gene)) %>% mutate(novelTranscript = ifelse(associated_transcript == "novel","Novel","Known"))

# recount support 
recountSupport <- fread(paste0(root_dir, "sqanti/sqanti_final_recount_support.txt.gz"))

# further filtering of 10 reads, 10 samples 
# more stringent filtering beyond 2 reads, 2 samples
# whole+targeted, nreads, nsamples
class.files$glob_targ_SQ <- class.files$glob_targ_SQ %>% filter(nreads >= 10 & nsamples >= 10)
class.files$glob_targ_SQ_default <- class.files$glob_targ_SQ_default %>% filter(nreads >= 10 & nsamples >= 10)
# glob_SQ: whole_nreads, whole_nsamples
class.files$glob_SQ_default <- class.files$glob_targ_SQ_default %>% filter(whole_nreads >= 10 & whole_nsamples >= 10)
class.files$glob_SQ_default_annoGene <- class.files$glob_targ_SQ_default_annoGene %>% filter(whole_nreads >= 10 & whole_nsamples >= 10) 
class.files$glob_SQ_annoGene <- class.files$glob_SQ_annoGene %>% filter(whole_nreads >= 10 & whole_nsamples >= 10) 
class.files$glob_SQ <- class.files$glob_SQ %>% filter(whole_nreads >= 10 & whole_nsamples >= 10)
# targeted: targeted_nreads, targeted_nsamples
class.files$targ_SQ <- class.files$glob_SQ_annoGene %>% filter(targeted_nreads >= 10 & targeted_nsamples >= 10)

# only keep if full support
recountSupport <- recountSupport %>% mutate(ToKeep = ifelse(propsupportedjunctions == 1, TRUE, FALSE))
# To keep ISM and FSM regardless of recount support
recountSupport <- recountSupport %>% mutate(ToKeep = ifelse(structural_category %in% c("FSM","ISM"), TRUE, ToKeep))
# unique(recountSupport[recountSupport$structural_category == "FSM","ToKeep"])
# unique(recountSupport[recountSupport$structural_category == "ISM","ToKeep"])
# unique(recountSupport[!recountSupport$structural_category %in% c("FSM","ISM") & recountSupport$ToKeep == TRUE, "propsupportedjunctions"])
supportedRecountSupport <- recountSupport[recountSupport$ToKeep == TRUE, ]
class.files <- lapply(class.files, function(x) x %>% filter(isoform %in% supportedRecountSupport$isoform))

# check all transcripts in glob_SQ_annoGene in glob_SQ
if(length(setdiff(class.files$glob_SQ_annoGene$isoform, class.files$glob_SQ$isoform)) != 0){
  message("Error: isoforms in class.files$glob_SQ_annoGene not in class.files$glob_SQ")
}

annoGenesStats <- list(
  novelTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript == "novel",],
  annoTrans = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$associated_transcript != "novel",],
  NIC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NIC",],
  NNC = class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$structural_category == "NNC",]
)

## Protein level 
class.files$protein_filtered = fread(paste0(root_dir,"/proteomics/Whole.sqanti_protein_classification.tsv"),data.table = F)
# keep only 10 reads, 10 samples isoforms
class.files$protein_filtered <- class.files$protein_filtered %>% filter(pb %in% class.files$glob_targ_SQ$isoform)
# annotate ensembl to gene symbol
class.files$protein_filtered$associated_gene <- class.files$protein_filtered$pr_gene
class.files$protein_filtered <- convert_ensembl_to_genename(unique(class.files$protein_filtered))

# cpat
cpat <- fread(paste0(root_dir,"/proteomics/WholeTargeted_filtered_finalversion.ORF_prob.best.tsv"), stringsAsFactors = F, data.table = F)

## prenatal vs postnatal
class.files$glob_SQ_annoGene_prenatal <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$preReads >= 1,]
class.files$glob_SQ_annoGene_postnatal <- class.files$glob_SQ_annoGene[class.files$glob_SQ_annoGene$postReads >= 1,]

## transcripts annotated to known protein-coding genes
class.files$glob_SQ_proteinGenes <- class.files$glob_SQ[class.files$glob_SQ$associated_gene %in% ProteinCodingGenes,]

## intergenic and fusion (from K.Chundru with strict filters)
fusion <- read.csv(paste0(root_dir, "/sqanti/fusion.csv")) %>% filter(isoform %in% class.files$glob_SQ$isoform)
intergenic <- read.csv(paste0(root_dir, "/sqanti/intergenic.csv")) %>% filter(isoform %in% class.files$glob_SQ$isoform)

## generate FL reads for QC report
wholeTargetedDemux <- class.files$glob_targ_SQ %>% select(isoform, contains("Whole"), -whole_nreads, -whole_nsamples) 
colnames(wholeTargetedDemux)[1] <- "id"
#write.csv(wholeTargetedDemux, paste0(output_dir, "/wholeTargeted_demux_fl_count.csv"), row.names = F)

# write output
#write.table(class.files$glob_SQ, paste0(root_dir, "sqanti/sqantifiltered_monoexonicfiltered_10reads10samplesRecount_whole_classification.txt"), sep = "\t", row.names = F, quote = F)
#write.table(class.files$glob_targ_SQ, paste0(root_dir, "sqanti/sqantifiltered_monoexonicfiltered_10reads10samplesRecount_classification.txt"), sep = "\t", row.names = F, quote = F)

# genrate gtf for recount
write.table(class.files$glob_targ_SQ$isoform, paste0(root_dir, "sqanti/sqantifiltered_monoexonicfiltered_10reads10samplesRecount_isoform.txt"), 
            sep = "\t", col.names = F, row.names = F, quote = F)
# gtf 
#originalGtf=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf
#cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/
#grep -f sqantifiltered_monoexonicfiltered_10reads10samplesRecount_isoform.txt ${originalGtf} > sqantifiltered_monoexonicfiltered_10reads10samplesRecount.gtf

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


## -------------- differential transcript expression ----------------

WholeDTE <- list(
  sex = vroom(paste0(root_dir,"DTE/DESeq2_whole_transcript_sex_resSig.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(root_dir,"DTE/DESeq2_whole_transcript_development_resSig.csv"),delim = ",",show_col_types = FALSE)
)
WholeDTE <- lapply(WholeDTE, function(x) x %>% mutate(dirAcrossDev = ifelse(log2FoldChange < 0 , "upregulated", "downregulated")))
WholeDTE$sex <- merge(WholeDTE$sex, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)
WholeDTE$age <- merge(WholeDTE$age, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)

WholeDTEAll <- list(
  sex = vroom(paste0(root_dir,"DTE/DESeq2_whole_transcript_sex_resAll.csv"), delim = ",",show_col_types = FALSE)
)
WholeDTEAll$sex <- merge(WholeDTEAll$sex, class.files$glob_targ_SQ[,c("isoform","chrom")],by.x="isoform",all.x=T)


## -------------- differential gene expression -------------

WholeDESeqGeneSig <- list(
  sex = as.data.frame(fread(paste0(root_dir,"DGE/DESeq2_whole_gene_sex_resSig.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"DGE/DESeq2_whole_gene_development_resSig.csv")))
)

WholeDESeqGene <- list(
  sex = as.data.frame(fread(paste0(root_dir,"DGE/DESeq2_whole_gene_sex_resAll.csv"))),
  age = as.data.frame(fread(paste0(root_dir,"DGE/DESeq2_whole_gene_development_resAll.csv")))
)

## -------------- differential isoform usage ----------------

load(file = paste0(root_dir,"DIU","/DIUSig.RData"))
load(file = paste0(root_dir,"DIU","/DIU.RData"))
DIUSig$wholeAge <- DIUSig$wholeAge %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
                                                 DGE_Sex = ifelse(Gene %in% WholeDESeqGeneSig$sex$associated_gene,TRUE,FALSE))
DIUSig$wholeSex <- DIUSig$wholeSex %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
                                              DGE_Sex = ifelse(Gene %in% WholeDESeqGeneSig$sex$associated_gene,TRUE,FALSE))
DIUSig <- lapply(DIUSig, function(x) x[complete.cases(x),])
DIUnormExp = list(
  GNAS_sex = paste0(root_dir, "DIU/GNAS_normalised_expression.txt"),
  GMP6A_group = paste0(root_dir, "DIU/GPM6A_normalised_expression.txt"),
  MORF4L2_group = paste0(root_dir, "DIU/MORF4L2_normalised_expression.txt"),
  HSP90AA1_sex = paste0(root_dir, "DIU/HSP90AA1_normalised_expression.txt")
)

## -------------- normalized counts ----------------

# differentially expressed transcripts
Exp <- list(
  whole_group = vroom(paste0(root_dir,"DTE/DESeq2_whole_development_normSig.csv"),delim = ",", show_col_types = FALSE),
  whole_sex = vroom(paste0(root_dir,"DTE/DESeq2_whole_sex_normSig.csv"),delim = ",", show_col_types = FALSE)
)
Exp <- lapply(Exp, function(x) merge(x, phenotype, by="sample"))

ExpGene <- list(
  whole_group = vroom(paste0(root_dir,"DGE/DESeq2_whole_development_normAll.csv"),delim = ",", show_col_types = FALSE),
  whole_sex = vroom(paste0(root_dir,"DGE/DESeq2_whole_sex_normAll.csv"),delim = ",", show_col_types = FALSE)
)

# normalised counts for all transcripts
normWhole <- fread(paste0(root_dir,"DTE/DESeq2_whole_development_normAll.csv"))
WholePreNorm <- normWhole %>% filter(sample %in% preWhole) %>% group_by(isoform) %>% tally(normalised_counts)
WholePostNorm <- normWhole %>% filter(sample %in% postWhole) %>% group_by(isoform) %>% tally(normalised_counts)
normWholeIsoform <- merge(WholePreNorm %>% `colnames<-`(c("isoform", "normPre")) , WholePostNorm %>% `colnames<-`(c("isoform", "normPost")), by = "isoform")

## -------------- gtf ----------------

gtf <- list()
gtf$ONT10.5139.1910 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/ONT10.5139.1910.gtf"))) %>% mutate(gene_id = "ADD3", transcript_id = "ONT10.5139.1910")
gtf$ONT18.5258.1932 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/ONT18.5258.1932.gtf"))) %>% mutate(gene_id = "MBP", transcript_id = "ONT18.5258.1932")
gtf$ONT2.10213.11813 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/ONT2.10213.11813.gtf"))) %>% mutate(gene_id = "CHN1", transcript_id = "ONT2.10213.11813")
gtf$GNAS <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/GNAS.gtf"))) 
gtf$MORfL2 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/MORf4L2.gtf"))) 
gtf$GPM6A <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/GPM6A.gtf"))) 
gtf$DLGAP5 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/DLGAP5.gtf"))) 
gtf$FOXP2 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/FOXP2.gtf"))) 
gtf$DAGLA <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/DAGLA.gtf"))) 
gtf$CACNA1G <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/CACNA1G.gtf"))) 
gtf$RPS27A <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/RPS27A.gtf"))) 
gtf$EIF2S3 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/EIF2S3.gtf"))) 
gtf$TCF4 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/TCF4.gtf"))) 
gtf$HERC1 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/HERC1.gtf"))) 
gtf$GRIA3 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/GRIA3.gtf")))
gtf$ZMYM2 <- as.data.frame(rtracklayer::import(paste0(root_dir,"sqanti/ZMYM2.gtf")))

gtf$ref <- data.table::fread(paste0(root_dir,"utils/refExons.gtf")) %>% dplyr::rename("gene_id" = "gene_name") %>% mutate(type = "exon")
gtf <- lapply(gtf, function(x) x[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

gtf$merged <- rbind(gtf$ONT18.5258.1932[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$ONT10.5139.1910[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$ONT2.10213.11813[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$GNAS[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$MORfL2[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$GPM6A[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$FOXP2[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$DAGLA[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$DLGAP5[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$CACNA1G[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$RPS27A[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$EIF2S3[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$TCF4[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$HERC1[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$GRIA3[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$ZMYM2[,c("seqnames","strand","start","end","type","transcript_id","gene_id")],
                    gtf$ref[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

GI <- c("GRIN2A","GRIA3","SEPTIN4","RTN4","MBP","RPS4Y1","XIST","ADD3","CNTNAP2","ANKRD12","VXN","PKM","MORF4L2","GNAS","RSP27A","CHN1","GPM6A",
        "TCF4","HERC1","ZMYM2")
RefIsoforms <- lapply(GI, function(x) unique(gtf$ref[gtf$ref$gene_id == x & !is.na(gtf$ref$transcript_id), "transcript_id"]))
names(RefIsoforms) <- GI


## -------------- long read proteoogenomics ----------------

# created from 7_process_proteomics_input.R
load(paste0(root_dir,"proteomics/proteinInputWhole.RData"))

# filter protein input to only the class files (pb_accs as this was the original pb_accs before being collapsed)
proteinInput$t2p.collapse <- proteinInput$t2p.collapse %>% filter(pb_accs %in% class.files$glob_SQ$isoform)

proteinInput$filtered <- fread(paste0(root_dir,"proteomics/Whole.classification_filtered.tsv"), data.table = F)


## -------------- mass spectrometry ----------------
# All peptides
load(paste0(root_dir,"proteomics/AllPeptides.RData"))

# Novel peptides
load(paste0(root_dir,"proteomics/NovelPeptides.RData"))


## -------------- ERCC ----------------

ercc.class.file <- fread(paste0(root_dir,"ERCC/AlessiaERCC_collapsed_RulesFilter_result_classification.txt"), sep = "/t") %>% filter(filter_result == "Isoform")
ercc.demux <- fread(paste0(root_dir,"ERCC/demux_fl_count.csv"))
ercc.stats <- fread(paste0(root_dir,"ERCC/ERCC_Stats.csv"))

## -------------- comparison with Leung et al.(2021) dataset ----------------

# gffcompare output 
gfftmapComparisons <- list(
  cellReports = data.table::fread(paste0(root_dir, "overlapDatasets/sfari_PacBio.HumanCTX.collapsed_classification.filtered_lite.gtf.tmap"), data.table = FALSE),
  directRNA = data.table::fread(paste0(root_dir, "overlapDatasets/sfari_dRNA.sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf.tmap"), data.table = FALSE)#,
  #BDRNatureComms = data.table::fread(paste0(root_dir, "overlapDatasets/sfari_BDR.ontBDR_collapsed.filtered_counts_filtered.gtf.tmap"), data.table = FALSE),
  #Patowary = data.table::fread(paste0(root_dir, "overlapDatasets/sfari_Patowary.cp_vz_0.75_min_7_recovery_talon_corrected.gtf.tmap"), data.table = FALSE),
  #PatowaryHerberle = data.table::fread(paste0(root_dir,"overlapDatasets/Heberle_Patowary.cp_vz_0.75_min_7_recovery_talon_corrected.gtf.tmap"), data.table = FALSE),
  #Herberle = data.table::fread(paste0(root_dir,"overlapDatasets/Bamford_Heberle.extended_annotations.gtf.tmap"), data.table = FALSE),
  #HerberleBamford = data.table::fread(paste0(root_dir,"overlapDatasets/Heberle_Bamford.sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf.tmap"), data.table = FALSE)
)

# from zenodo
humanCTX <- read.table(paste0(root_dir, "overlapDatasets/HumanCTX.collapsed_classification.filtered_lite_classification.txt"), header = TRUE)
humanCTX$totalFL <- humanCTX %>% dplyr::select(contains("FL.")) %>% apply(.,1,sum)
# copied from here: /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/9_sqanti_final
directRNA <- fread(paste0(root_dir, "overlapDatasets/sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt"), data.table = FALSE)
# copied from here: /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3
BDRNatureComms <- data.table::fread(paste0(root_dir, "overlapDatasets/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt"),data.table = FALSE)
PatowaryCTX <- data.table::fread(paste0(root_dir, "overlapDatasets/cp_vz_0.75_min_7_recovery_talon_classification.txt"), data.table = FALSE)
HerberleCTX <- data.table::fread(paste0(root_dir, "overlapDatasets/fullLengthCounts_transcript.txt"), data.table = FALSE) %>% mutate(isoform = TXNAME)

# subset global targeted dataset to 20 AD target genes
AD20TargetGenes <- c("ABCA1", "PICALM", "SORL1", "FUS", "MAPT", "RHBDF2", "ABCA7", "APOE", "CD33", "BIN1", "TARDBP", "APP", "SNCA", "TREM2", "FYN",
  "VGF", "PTK2B","CLU", "ANK1", "TRPA1")
class.files$glob_targ_SQ_20AD <- class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene %in% AD20TargetGenes,]


## -------------- FICLE ----------------

# read in FICLE transcript_classifications.csv merged
finalTranscriptClassification <- distinct(fread(paste0(root_dir,"ficle/all_final_transcript_classifications.csv"), data.table = F))
finalTranscriptClassification <- finalTranscriptClassification %>% filter(isoform %in% class.files$glob_SQ$isoform)
finalTranscriptClassificationTranscript <- merge(finalTranscriptClassification, class.files$glob_SQ[,c("isoform","associated_gene","preReads", "postReads","DevStatus")], by = "isoform")
finalTranscriptClassificationTranscript <- merge(finalTranscriptClassificationTranscript, normWholeIsoform, by = "isoform", all.x = T)

# numeric for columns responding to AS events
finalTranscriptClassificationTranscript$isoform <- as.factor(finalTranscriptClassificationTranscript$isoform)
finalTranscriptClassificationTranscript$associated_gene <- as.factor(finalTranscriptClassificationTranscript$associated_gene)
finalTranscriptClassificationTranscript$DevStatus <- as.factor(finalTranscriptClassificationTranscript$DevStatus)
finalTranscriptClassificationTranscript <- finalTranscriptClassificationTranscript %>% mutate_if(is.character, as.numeric)
#length(unique(finalTranscriptClassificationGene$associated_gene))

# aggregate the number of splicing events per column 
finalTranscriptClassificationGene <- aggregate(. ~ associated_gene, finalTranscriptClassificationTranscript %>% select(-isoform, - DevStatus), sum)

finalTranscriptClassificationGeneDev <- aggregate(. ~ associated_gene + DevStatus, 
                                               finalTranscriptClassificationTranscript %>% select(-isoform), sum)
finalTranscriptClassificationGeneDev <- reshape2::melt(finalTranscriptClassificationGeneDev, 
                                                     variable.name = "AS", value.name = "Frequency", id = c("associated_gene","DevStatus"))

## ----- selected disease genes ------
monoAllelicDDP = c('ADNP','ANK2','ANKRD11','ARID1B','ASH1L','AUTS2','BCL11A','BCL11B','CACNA1C','CACNA1G','CACNA1H','CDH1','CDH2','CDK13','CHD2','CHD8','CLCN3','CNOT1','DLG4','DMD','DYRK1A','EIF2S3','EP300','FMR1','FOXP1','FOXP2','GABBR2','GATA3','GATA4','GATA6','GLMN','GRIA3','GRIK2','GRIN1','GRIN2A','GRIN2B','H1-4','HNF1B','HNF4A','IL1RAPL1','ITPR1','KAT6B','KDM6A','KDM6B','KIRREL3','KMT2D','LMNA','MAGI2','MECP2','MED13L','MEIS2','MNX1','NF1','NFIA','NFIB','NLGN3','NLGN4X','NRXN1','OPA1','PAX6','PBX1','POGZ','POLA1','PTEN','QRICH1','RANBP2','RBFOX1','RPL10','SCN1A','SCN2A','SCN8A','SET','SETD1A','SHANK3','SHH','SLC9A9','SMARCE1','SON','SOX9','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','UBE3A','USP7','WFS1','ZBTB20','ZEB2','ZMYM2')

biallelicDDP = c('ASPM','CACNA1G','CHL1','CISD2','CLCN3','CNTNAP2','CTNNA2','DCC','EIF2AK3','EOMES','GLIS3','GPC6','GRIK2','GRIN1','GRIN2A','GRM1','HADH','HERC1','HPSE2','ITPR1','KDM6B','KIAA1109','LMNA','LRBA','MCPH1','NRXN1','NTNG1','ONECUT1','PDIA6','PDX1','PEX16','PPP1R15B','PTF1A','QARS1','RELN','RFX6','SLC39A8','SLF2','TRMT10A','TSPEAR','UBE3B','VPS13B','WFS1','WWOX','ZBTB16','ZFP57')

GWAS = c("ACTR1B", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CSMD1", "CTD-2008L17.2", "CUL9", "DCC", "DLGAP2", "DPYD", "EMX1", "ENOX1", "EPN2", "EYS", "FURIN", "GABBR2", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "IL1RAPL1", "IMMP2L", "IRF3", "KIAA1549", "KLF6", "LINC00320", "LINC01088", "LRRC4B", "MAD1L1", "MAN2A1", "MAPT", "MSI2", "NAB2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCGF3", "PCNXL3", "PDE4B", "PJA1", "PLCH2", "PTPRD", "R3HDM2", "RP11-399D6.2", "RP11-507B12.2", "SGCD", "SLC39A8", "SLC4A10", "SNAP91", "SP4", "THAP8", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WSCD2", "ZNF804A", "ZNF823", "ZNF835")

SFARI = c('ADNP','AGAP1','ANK2','ANKRD11','ARID1B','ASH1L','ASPM','AUTS2','BCL11A','CACNA1C','CACNA1G','CACNA1H','CADPS','CADPS2','CD38','CDH13','CDH2','CDK13','CELF6','CHD2','CHD8','CLTCL1','CNOT1','CNTN4','CNTN6','CNTNAP2','CSMD1','CTNNA2','CYFIP1','DAGLA','DCC','DISC1','DLG4','DLGAP2','DMD','DPYD','DRD2','DRD3','DYRK1A','ELP4','EP300','FBXO40','FMR1','FOXP1','FOXP2','GABBR2','GPC6','GRIA3','GRIK2','GRIK3','GRIN1','GRIN2A','GRIN2B','H1-4','HERC1','IL1RAPL1','IMMP2L','ITPR1','KAT6B','KATNAL2','KCTD13','KDM6A','KDM6B','KIRREL3','LRBA','MAPT','MCM4','MCPH1','MECP2','MED13L','MEIS2','MET','NEGR1','NF1','NFIA','NFIB','NKX2-2','NLGN1','NLGN2','NLGN3','NLGN4X','NR3C2','NRXN1','NTNG1','NXPH1','OXTR','PARD3B','PAX6','PBX1','PCDH10','PCDH9','PHB','PJA1','POGZ','PRKN','PTEN','PTPRT','QRICH1','RBFOX1','RELN','RPL10','RUNX1T1','SCN1A','SCN2A','SCN8A','SET','SETD1A','SEZ6L2','SHANK3','SLC4A10','SLC9A9','SON','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','TSHZ3','UBE3A','USP7','VPS13B','WWOX','ZBTB16','ZBTB20','ZMYM2','ZNF18','ZNF804A')

schemaGenes <- read.table(paste0(root_dir, "metadata/schema_genes.txt"), col.names = c("Gene"))

selectedTargetGenes <- unique(c(as.character(monoAllelicDDP), as.character(SFARI), as.character(schemaGenes$Gene), as.character(biallelicDDP), GWAS))
length(selectedTargetGenes)
setdiff(selectedTargetGenes, TargetGene)
setdiff(selectedTargetGenes, TargetGeneS2$V1)
setdiff(TargetGeneS2$V1, selectedTargetGenes)
