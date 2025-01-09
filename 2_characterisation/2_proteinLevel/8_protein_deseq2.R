## ---------- Script -----------------
##
## Purpose: perform differential analysis on peptides predicted from CPAT
## Transcript level separate analysis
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
##
## ---------- Notes ------------------


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))


## ---------- source functions -----------------

LOGEN <- "C:/Users/sl693/Dropbox/Scripts/LOGen/"
source(paste0(LOGEN, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN, "differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))


## ---------- input -----------------

# directory names
root_dir <- "C:/Users/sl693/Dropbox/Scripts/SFARI_developmentalgenomics/data/"
output_dir <- "C:/Users/sl693/Dropbox/Scripts/SFARI_developmentalgenomics/output/"

# classfiles
load(file = paste0(root_dir,"sqanti/sqantifiltered_monoexonicfiltered_2reads2samples.RData"))

# proteomics input 
load(file = paste0(root_dir,"proteomics/proteinInputWhole.RData"))

input <- list()
# phenotype
input$wholePhenotype <- fread(paste0(root_dir, "metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% 
  mutate(group = factor(group, levels = c("Prenatal","Postnatal")), col = paste0(sample,"_",group)) %>% 
  filter(grepl("Whole",sample))
input$classfiles <- class.files$glob_SQ %>% select(isoform, associated_gene, associated_transcript, structural_category, subcategory, 
                                                   contains("Whole", ignore.case = FALSE))
row.names(input$classfiles) <- input$classfiles$isoform


# phenotype
phenoFiles <- list(
  wholeGroup = input$wholePhenotype,
  wholeSex = input$wholePhenotype %>% mutate(group = sex)
)

### protein 
# aggregate sum by same peptide sequence
counts <- class.files$glob_SQ %>% select(isoform, contains("Whole", ignore.case = FALSE))
pFL <- proteinInput$t.class.files %>% dplyr::select(base_acc, contains("Whole", ignore.case = FALSE)) 
pFLsum <- aggregate(. ~ base_acc, pFL, sum)

# datawrangle for input to run_DESeq2()
expressionFiles <- list(
  iso = input$classfiles %>% dplyr::select(starts_with("Whole")),
  isoProtein = pFLsum %>% tibble::column_to_rownames(., var = "base_acc")
)
row.names(expressionFiles$iso) <- input$classfiles$isoform


## ---------- ONT: Creating DESeq2 object and analysis -----------------

ResTran <- list(
  pWald = run_DESeq2(test="Wald",expressionFiles$isoProtein,phenoFiles$wholeGroup,threshold=10,controlname="Prenatal",design="case_control",groupvar="factor")
)

ResTran$pWald$annoRes <- merge(ResTran$pWald$res_Wald[ResTran$pWald$res_Wald$padj < 0.05, ],
      input$classfiles[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory")], by = "isoform")

## ---------- Output -----------------

saveRDS(ResTran, file = paste0(output_dir, "/IsoProtein_Whole_DESeq2TranscriptLevel.RDS"))
