#!/usr/bin/env Rscript
## ----------Script-----------------
## 
## SK.Leung: filter the sqanti files generated per chromosome by 
##           keeping only the transcripts kept from sqanti filtering
##           removing mono-exonic intergenic transcripts
##           removing mono-exonic transcripts within multi-exonic genes
## ---------------------------------

library("stringr")

## ------------ directory names --------------- 

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))


## -------------- Final classification files -------------

dirnames = list(
	sqanti = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax_merged/",
	utils = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_utils/",
  output ="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_output/"
)
class.names.files = list.files(path = dirnames$sqanti, pattern = "RulesFilter_result_classification.txt", full = T)

class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
names(class.files) <- list.files(path = dirnames$sqanti, pattern = "RulesFilter_result_classification.txt")
save(class.files, file = paste0(dirnames$output,"sqantifiltered_classification.RData"))

# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons != "Intergenic_1"))

# remove mono-exonic transcripts 
refExonNum <- read.csv(paste0(dirnames$utils,"gencode.v40.annotation.numExon.csv"))
exonicClass <- refExonNum %>% 
  group_by(Gene, GeneName, GeneType) %>% 
  dplyr::summarise(maxNumExon = max(NumExon)) %>% 
  mutate(monoExonic = ifelse(maxNumExon == 1, TRUE, FALSE))

refGeneType <- read.table(paste0(dirnames$utils,"gencode.v38.annotation.geneannotation.txt"), header = T)

# identify multi-exonic genes
multiExonicGenes <- unique(exonicClass[exonicClass$monoExonic == FALSE,][["GeneName"]])
# TRUE = monotranscript within multi-exonic gene
class.files <- lapply(class.files, function(x) x %>% mutate(monomulti = ifelse(exons == 1 & associated_gene %in% multiExonicGenes, TRUE,FALSE)))
# filter FALSE to retain multi-transcripts within multi-exonic gene, and mono-transcripts within mono-exonic gene
mono.multi.class.files <- lapply(class.files, function(x) x %>% filter(monomulti == TRUE))
class.files <- lapply(class.files, function(x) x %>% filter(monomulti == FALSE))
save(class.files, file = paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_classification.RData"))

merged.class.files <- bind_rows(class.files, .id = "filename")
merged.class.files <- merged.class.files %>% select(-filename)
write.table(merged.class.files, file = paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_classification.txt"), quote=F, sep = "\t", row.names= F)
write.table(merged.class.files$isoform, file = paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_ID.txt"), quote=F, sep = "\t", row.names= F, col.names = F)