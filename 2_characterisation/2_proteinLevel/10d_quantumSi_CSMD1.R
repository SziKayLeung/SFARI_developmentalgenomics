suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("vroom"))
suppressMessages(library("tidyr"))
suppressMessages(library("stringr"))


## ------------ directory names --------------- 

LOGEN <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/LOGen"
root_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/PaperZenodo/"
output_dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/SFARI/Output"
SC_ROOT= "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))

## ------------- Phenotype files -------------------

phenotype <- fread(paste0(root_dir, "metadata/WholeTargetedphenotype_fixedsex.csv"),data.table=F, stringsAsFactors=F) %>% mutate(time = age)
phenotype <- phenotype %>% mutate(type = ifelse(grepl("Targeted",sample),"Targeted","Whole"),
                                  sampleID = gsub("^Targeted", "", sample)) %>% mutate(sampleID = gsub("^Whole","", sampleID))
phenotype <- phenotype %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")), 
                                  col = paste0(sample,"_",group),
                                  sex = factor(sex, levels = c("M","F"))) 

# manifest 
manifest <- fread(paste0(root_dir, "metadata/WholeTargetedphenotype_manifest.csv"),data.table=F, stringsAsFactors=F)

## -------------- Final classification files ------------- 

load(file = paste0(root_dir,"sqanti/sqantifiltered_monoexonicfiltered_2reads2samples.RData"))
View(class.files$glob_targ_SQ[class.files$glob_targ_SQ$isoform %in% c("ONT8.85.5521","ONT8.85.5571"),
                         c("isoform","chrom","strand", "length", "exons", "structural_category", "associated_gene",
                           "all_canonical","within_polya_site","polyA_motif_found",
                           "whole_nreads","whole_nsamples","targeted_nreads","targeted_nsamples")])

View(class.files$glob_targ_SQ[class.files$glob_targ_SQ$isoform == "ONT3.117.12367",
                              c("isoform","chrom","strand", "length", "exons", "structural_category", "associated_gene",
                                "all_canonical","within_polya_site","polyA_motif_found",
                                "whole_nreads","whole_nsamples","targeted_nreads","targeted_nsamples")])


WholeDTE <- list(
  sex = vroom(paste0(root_dir,"DTE/DESeq2_whole_transcript_sex_resSig.csv"),delim = ",",show_col_types = FALSE),
  age = vroom(paste0(root_dir,"DTE/DESeq2_whole_transcript_development_resSig.csv"),delim = ",",show_col_types = FALSE)
)
WholeDTE <- lapply(WholeDTE, function(x) x %>% mutate(dirAcrossDev = ifelse(log2FoldChange < 0 , "upregulated", "downregulated")))
Exp <- list(
  whole_group = vroom(paste0(root_dir,"DTE/DESeq2_whole_development_normSig.csv"),delim = ",", show_col_types = FALSE),
  whole_sex = vroom(paste0(root_dir,"DTE/DESeq2_whole_sex_normSig.csv"),delim = ",", show_col_types = FALSE)
)
Exp <- lapply(Exp, function(x) merge(x, phenotype, by="sample"))


plot_trans_exp_individual("ONT8.85.5521",class.files$glob_SQ,Exp$whole_group,"group", 
                          sqrt=TRUE, colourdots = alpha("#00BFC4",0.3))
ONT8.85.5521Exp <- Exp$whole_group[Exp$whole_group$isoform == "ONT8.85.5521",c("sample","normalised_counts","group","sex","age","RIN","time")] %>% 
  filter(normalised_counts > 0) %>% arrange(normalised_counts)
View(class.files$glob_targ_SQ[class.files$glob_targ_SQ$isoform == "ONT8.85.5521",] %>% select(isoform, contains("Whole")) %>% 
  select(-whole_nsamples, -whole_nreads) %>% 
  reshape2::melt(variable.name = "sample", value.name = "FLreads") %>% filter(FLreads > 0) %>% 
  merge(., ONT8.85.5521Exp, by = "sample") %>% 
  arrange(-FLreads))

plot_trans_exp_individual("ONT3.117.12367",class.files$glob_SQ,Exp$whole_group,"group", 
                          sqrt=TRUE, colourdots = alpha("#00BFC4",0.3))
ONT3.117.12367Exp <- Exp$whole_group[Exp$whole_group$isoform == "ONT3.117.12367",c("sample","normalised_counts","group","sex","age","RIN","time")] %>% 
  filter(normalised_counts > 0) %>% arrange(normalised_counts)
View(class.files$glob_targ_SQ[class.files$glob_targ_SQ$isoform == "ONT3.117.12367",] %>% select(isoform, contains("Whole")) %>% 
       select(-whole_nsamples, -whole_nreads) %>% 
       reshape2::melt(variable.name = "sample", value.name = "FLreads") %>% filter(FLreads > 0) %>% 
       merge(., ONT3.117.12367Exp, by = "sample") %>% 
       arrange(-FLreads))


load(paste0(root_dir,"proteomics/proteinInputWhole.RData"))
proteinInput$t.class.files[proteinInput$t.class.files$isoform == "ONT3.117.12367",]


## -------------- plot isoform usage ------------- 

# take only the samples of interest and calculate the normalised counts (TPM)
norm.class.files <- class.files$glob_targ_SQ[, c("isoform","associated_gene","WholeS1401","WholeDH3671","Whole45988")] %>% 
  mutate(across(where(is.numeric), ~ . / sum(.) * 1000000, .names = "norm_{.col}")) %>%
  select(isoform,associated_gene,contains("norm"))

# subset to CMSD1 only and calculate percentage (isoform fraction)
CSMD1 <- norm.class.files[norm.class.files$associated_gene == "CSMD1",]
CSMD1 <- CSMD1 %>%  mutate(across(where(is.numeric), ~ . / sum(.) * 100, .names = "perc_{.col}")) %>% select(isoform,contains("perc"))
CSMD1_melt <- reshape2::melt(CSMD1, value.name = "perc", variable.name = "sample")  

# minor isoforms = < 1%
# sum percentage across the minor isoforms
minorLessThan1 <- CSMD1_melt %>% filter(perc < 1) 
minor <- CSMD1 %>% filter(isoform %in% minorLessThan1$isoform) %>% select(isoform, contains("perc"))
minorgrouped <- minor %>% tibble::column_to_rownames(., var = "isoform") %>% apply(., 2, sum) %>% as.data.frame()
colnames(minorgrouped) <- c("minor")

# subset the major isoforms
major <- CSMD1_melt %>% filter(!isoform %in% minorLessThan1$isoform) 
minorgrouped <- minorgrouped %>% tibble::rownames_to_column(., var = "sample") %>% 
  reshape2::melt(., value.name = "perc", variable.name = "isoform")

# merge 
df <- rbind(major,minorgrouped)
df <- df %>% mutate(sample = word(sample, c(3), sep = fixed("_")))

ggplot(df, aes(x = sample, y = perc, fill = isoform)) + geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "Proportion", x = "Sample")

