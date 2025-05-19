## ---------- Script -----------------
##
## Purpose: identify differentially expressed novel transcripts (translated to novel protein isoforms) between prenatal and postnatal datasets
## Based from results in 8_protein_deseq2.R
## Top 2 differentially expressed peptides between prenatal and postnatal 
## Top 20 differentially expressed peptides associated with SFARI genes
##
## ---------- Notes ------------------
## identify peptides for Quantum-Si predicted from long-read sequencing sfari datasets
## n = 5,384 novel transcripts/peptides
## 
## differential expression (8_protein_deseq2.R) was performed using base_acc from original long-read proteogenomics pipeline (not the representative most abundant transcript)
## so can directly filter using the same base_acc from the sqanti classificaiton file generated from long-read proteogenomics pipeline

## ---------- packages -----------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

## ---------- input -----------------

# directory names
root_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/SFARI_developmentalgenomics/data/"
output_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/PaperZenodo/"

# differentially expressed transcripts collapsed (output of 8_protein_deseq2.R)
ResTran <- readRDS(paste0(output_dir, "/IsoProtein_Whole_DESeq2TranscriptLevel.RDS"))

# classification file from running long-read proteo-genomics pipeline
protein.class.files <- fread(paste0(root_dir, "proteomics/Whole.sqanti_protein_classification.tsv"))

# protein exploration (7_process_proteogenomics_exploration.R)
load(file = paste0(root_dir,"proteomics/proteinInputWhole.RData"))

## ----------------------------------

# 1. identify top-ranked differentially expressed novel transcripts (NIC, NNC only)
ResTran$pWald$annoResNovel <- ResTran$pWald$annoRes %>% arrange(padj) %>% filter(structural_category %in% c("NIC","NNC")) 

# 2. from those differentially expressed novel transcripts, find the ones that are translated to novel protein isoforms
# Note: novel transcript can translate to a known peptide sequence, novel transcript can be varying at the UTR
novelDifferentialProtein <- protein.class.files[protein.class.files$pb %in% ResTran$pWald$annoResNovel$isoform,]
novelDifferentialProtein <- novelDifferentialProtein[novelDifferentialProtein$pr_splice_cat %in% c("novel_in_catalog","novel_not_in_catalog"),]

# 3. filter for transripts that are not predicted for nonsense mediated decay, as unlikely to be detected in Quantum-Si
novelDifferentialProtein <- novelDifferentialProtein[novelDifferentialProtein$is_nmd == "FALSE",]

# 4. Plot transcripts differentially expressed between prenatal and postnatal
#ResTran$pWald$norm_counts[ResTran$pWald$norm_counts$isoform == "ONT2.10213.11813",] %>% 
#  merge(., phenoFiles$wholeGroup) %>% 
#  ggplot(., aes(x = group, y = normalised_counts)) + geom_point()


message("Number of novel differentially expressed transcripts predicted to code for novel protein: ", 
        nrow(ResTran$pWald$annoResNovel[ResTran$pWald$annoResNovel$isoform %in% novelDifferentialProtein$pb,]))

## ----------------------------------

## Top 20 novel differentially expressed transcripts associated to SFARI genes

# subset to SFARI gene only 
SFARI = c('ADNP','AGAP1','ANK2','ANKRD11','ARID1B','ASH1L','ASPM','AUTS2','BCL11A','CACNA1C','CACNA1G','CACNA1H','CADPS','CADPS2','CD38','CDH13','CDH2','CDK13','CELF6','CHD2','CHD8','CLTCL1','CNOT1','CNTN4','CNTN6','CNTNAP2','CSMD1','CTNNA2','CYFIP1','DAGLA','DCC','DISC1','DLG4','DLGAP2','DMD','DPYD','DRD2','DRD3','DYRK1A','ELP4','EP300','FBXO40','FMR1','FOXP1','FOXP2','GABBR2','GPC6','GRIA3','GRIK2','GRIK3','GRIN1','GRIN2A','GRIN2B','H1-4','HERC1','IL1RAPL1','IMMP2L','ITPR1','KAT6B','KATNAL2','KCTD13','KDM6A','KDM6B','KIRREL3','LRBA','MAPT','MCM4','MCPH1','MECP2','MED13L','MEIS2','MET','NEGR1','NF1','NFIA','NFIB','NKX2-2','NLGN1','NLGN2','NLGN3','NLGN4X','NR3C2','NRXN1','NTNG1','NXPH1','OXTR','PARD3B','PAX6','PBX1','PCDH10','PCDH9','PHB','PJA1','POGZ','PRKN','PTEN','PTPRT','QRICH1','RBFOX1','RELN','RPL10','RUNX1T1','SCN1A','SCN2A','SCN8A','SET','SETD1A','SEZ6L2','SHANK3','SLC4A10','SLC9A9','SON','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','TSHZ3','UBE3A','USP7','VPS13B','WWOX','ZBTB16','ZBTB20','ZMYM2','ZNF18','ZNF804A')

# list of IDs for subsequent downstream grep
TopRankedSFARINovelProteins <- ResTran$pWald$annoResNovel[ResTran$pWald$annoResNovel$isoform %in% novelDifferentialProtein$pb,] %>% 
  filter(associated_gene %in% SFARI) %>% 
  .[1:20,"isoform"]