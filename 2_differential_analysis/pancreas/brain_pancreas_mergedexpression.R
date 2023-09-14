## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
## 
## Prequisite: Run comparison of Brain targeted vs Pancreas targeted dataset (1_characterisation/3_brainvspancreas_characterisation.sh)
## Brain Fetal vs Pancreas Fetal: Targeted vs Targeted dataset 
##
##
##

## ---------- Packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))


## ---------- Config files -----------------

SOURCE_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/"
source(paste0(SOURCE_DIR,"sfari_differential.config.R"))


## ---------- Input files -----------------

Input_files <- list(
  gffcomp = gffcompare$brainvspancreas,
  Brain_exp = expression.files$brain_targeted,
  Pancreas_exp = expression.files$pancreas_targeted,
  Brain_gtf = gtf$brain_targeted,
  Pancreas_gtf = gtf$pancreas_targeted,
  Brain_junc = junction.files.names$brain_targeted, 
  Pancreas_junc = junction.files.names$pancreas_targeted
)


Read_files <- list(
  gffcomp = read.table(Input_files$gffcomp, header = T),
  Brain_exp = Input_files$Brain_exp %>% rename_all(paste0, "_Brain") %>% tibble::rownames_to_column(var = "ref"),
  Pancreas_exp = Input_files$Pancreas_exp %>% rename_all(paste0, "_Pancreas") %>% tibble::rownames_to_column(var = "query"),
  Brain_gtf = read.table(Input_files$Brain_gtf, sep = "\t") %>% mutate(isoform = word(word(V9,c(2),sep = fixed(" ")),c(1),sep = ";")) ,
  Pancreas_gtf = read.table(Input_files$Pancreas_gtf, sep = "\t") %>% mutate(isoform = word(word(V9,c(2),sep = fixed(" ")),c(1),sep = ";")),
  Brain_junc = read.table(Input_files$Brain_junc, header = T), 
  Pancreas_junc = read.table(Input_files$Pancreas_junc, header = T)
)

Read_files$gffcomp = read.table(Input_files$gffcomp, header = T)

## ---------- Merge expression -----------------

# common isoforms
gffcomp_match <- Read_files$gffcomp %>% filter(class_code == "=") %>% mutate(Brain_pancreas_ID = paste0(ref_id,"_",qry_id))

common_brain_exp <- Read_files$Brain_exp %>% filter(ref %in% gffcomp_match$ref_id) %>% 
  merge(., gffcomp_match[,c("ref_id","Brain_pancreas_ID")], by.x = "ref", by.y = "ref_id", all = T) %>% 
  tibble::column_to_rownames("Brain_pancreas_ID") %>% select(-ref) 
  
common_pancreas_exp <- Read_files$Pancreas_exp %>% filter(query %in% gffcomp_match$qry_id) %>% 
  merge(., gffcomp_match[,c("qry_id","Brain_pancreas_ID")], by.x = "query", by.y = "qry_id", all = T) %>% 
  tibble::column_to_rownames("Brain_pancreas_ID") %>% select(-query) 

if(nrow(common_brain_exp) == nrow(common_pancreas_exp)){
  common_merged_exp <- merge(common_brain_exp, common_pancreas_exp, by = 0)
}else{
  print("ERROR: Mismatched number of isoforms")
}


# unique isoforms 
unique_brain_exp <- Read_files$Brain_exp %>% filter(!ref %in% gffcomp_match$ref_id) %>% dplyr::rename("query" = "ref")
unique_pancreas_exp <- Read_files$Pancreas_exp %>% filter(!query %in% gffcomp_match$qry_id)

# if not commonly known isoforms (i.e. ENS...), then merge the unique isoforms into one dataset
if(length(grep("ENS", intersect(unique_pancreas_exp$query, unique_brain_exp$query))) == 0){
  # note need to create a query with dataset, as TALON reuses same TALON IDs for novel isoforms even if not the same isoform across datasets
  # query is therefore to generate unique identifier when merging
  unique_brain_exp = unique_brain_exp %>% mutate(query = paste0(query,"_","brain"))
  unique_pancreas_exp = unique_pancreas_exp %>% mutate(query = paste0(query,"_","pancreas"))
  unique_merged_exp = merge(unique_brain_exp, unique_pancreas_exp, by = "query", all = T)
}

# replace all NAs across expression with 0 
unique_merged_exp[is.na(unique_merged_exp)] <- 0


# merge common and unique isoforms into one dataset 
common_merged_exp <- common_merged_exp %>% dplyr::rename("isoform" = "Row.names")
unique_merged_exp <- unique_merged_exp %>% dplyr::rename("isoform" = "query")

final_exp <- rbind(common_merged_exp, unique_merged_exp)
nrow(common_merged_exp) + nrow(unique_merged_exp) == nrow(final_exp)


# keep only the fetal brain samples (# CONTROL)
fetal_samples <- metadata.files$brain_targeted[metadata.files$brain_targeted$Condition == "Control","SQ_ID"]
adult_samples <- metadata.files$brain_targeted[metadata.files$brain_targeted$Condition == "Case","SQ_ID"]
final_exp_fetal <- final_exp %>% select(-paste0(adult_samples,"_Brain"))
ncol(final_exp_fetal) == ncol(final_exp) - length(adult_samples)
setdiff(phenotype$ont_tissue$sample,colnames(final_exp_fetal))


# generate a phenotype file for tappAS (brain fetal targeted vs pancreas fetal targeted)
# keep duplicated samples
brain = data.frame(sample = fetal_samples, group = "CONTROL") 
pancreas = data.frame(sample = metadata.files$pancreas_targeted$SQ_ID, group = "CASE")
brain_pancreas_phenotype <- rbind(brain, pancreas) %>% mutate(sample = ifelse(group == "CONTROL", paste0(sample,"_Brain"), paste0(sample,"_Pancreas")))
length(brain_pancreas_phenotype[brain_pancreas_phenotype=="CONTROL"])
length(brain_pancreas_phenotype[brain_pancreas_phenotype=="CASE"])
cat("Total samples for expression analysis:", nrow(brain_pancreas_phenotype))
nrow(brain_pancreas_phenotype) + 1 == ncol(final_exp_fetal)

# merged classification file 
class.files$brain_pancreas <- list(
  brain_unique = class.files$brain_targeted %>% filter(!isoform %in% gffcomp_match$ref_id) %>% mutate(isoform = paste0(isoform,"_brain")), 
  pancreas_unique = class.files$pancreas_targeted %>% filter(!isoform %in% gffcomp_match$qry_id) %>% mutate(isoform = paste0(isoform,"_pancreas")),
  brain_pancreas = class.files$pancreas_targeted %>% filter(isoform %in% gffcomp_match$qry_id) %>% 
    merge(., gffcomp_match[,c("qry_id","Brain_pancreas_ID")], by.x = "isoform", by.y = "qry_id") %>% 
    mutate(isoform = Brain_pancreas_ID) %>% select(-Brain_pancreas_ID)
)

# merge
class.files$merged <- rbind(class.files$brain_pancreas$brain_unique,
      class.files$brain_pancreas$pancreas_unique,
      class.files$brain_pancreas$brain_pancreas)

# check that all the isoforms are kept and consistent across datasets
nrow(class.files$brain_pancreas$brain_unique) == nrow(unique_brain_exp)
nrow(class.files$brain_pancreas$pancreas_unique) == nrow(unique_pancreas_exp)
nrow(class.files$brain_pancreas$brain_pancreas) == nrow(common_merged_exp)


# subset gtf files 
Read_files$Brain_gtf <- Read_files$Brain_gtf %>% 
  mutate(isoform = word(word(V9,c(2),sep = fixed(" ")),c(1),sep = ";"),
         gene_id = paste0("gene_id \"", word(word(V9, c(2), sep = fixed(";")),c(3), sep = fixed(" ")), "\"",";")) 

Read_files$Pancreas_gtf <- Read_files$Pancreas_gtf %>% 
  mutate(isoform = word(word(V9,c(2),sep = fixed(" ")),c(1),sep = ";"),
         gene_id = paste0("gene_id \"", word(word(V9, c(2), sep = fixed(";")),c(3), sep = fixed(" ")), "\"",";")) 

subset_gtf <- list(
  brain_common = Read_files$Brain_gtf %>% filter(isoform %in% gffcomp_match$ref_id), 
  pancreas_common = Read_files$Pancreas_gtf %>% filter(isoform %in% gffcomp_match$qry_id) ,
  brain_unique = Read_files$Brain_gtf %>% filter(!isoform %in% gffcomp_match$ref_id) %>% 
    mutate(V10 = paste0("transcript_id \"", isoform, "_brain\"; ", gene_id)),
  pancreas_unique = Read_files$Pancreas_gtf %>% filter(!isoform %in% gffcomp_match$qry_id) %>% 
    mutate(V10 = paste0("transcript_id \"", isoform, "_pancreas\"; ", gene_id))
)

subset_gtf$merged_common <- merge(subset_gtf$pancreas_common, gffcomp_match[,c("qry_id", "Brain_pancreas_ID")], by.x = "isoform", by.y = "qry_id") %>% 
  mutate(V10 = paste0("transcript_id \"", Brain_pancreas_ID, "\"; ", gene_id))

cols = paste0("V",c(seq(1,8),10))

merged_gtf <- rbind(subset_gtf$brain_unique[cols], subset_gtf$pancreas_unique[cols], subset_gtf$merged_common[cols])
length(unique(merged_gtf$V10)) == length(unique(subset_gtf$brain_unique$V10)) + length(unique(subset_gtf$pancreas_unique$V10)) + length(unique(subset_gtf$merged_common$V10))

# check consistency
nrow(final_exp) == length(unique(merged_gtf$V10)) 
length(unique(merged_gtf$V10)) == nrow(class.files$merged)


# merged junction file 
junc.files <- list(
  brain_unique = Read_files$Brain_junc %>% filter(!isoform %in% gffcomp_match$ref_id) %>% mutate(isoform = paste0(isoform,"_brain")), 
  pancreas_unique = Read_files$Pancreas_junc %>% filter(!isoform %in% gffcomp_match$qry_id) %>% mutate(isoform = paste0(isoform,"_pancreas")),
  brain_pancreas = Read_files$Pancreas_junc %>% filter(isoform %in% gffcomp_match$qry_id) %>% 
    merge(., gffcomp_match[,c("qry_id","Brain_pancreas_ID")], by.x = "isoform", by.y = "qry_id") %>% 
    mutate(isoform = Brain_pancreas_ID) %>% select(-Brain_pancreas_ID)
)


# merge
junc.files$merged <- rbind(junc.files$brain_unique,
                           junc.files$pancreas_unique,
                           junc.files$brain_pancreas)

setdiff(unique(junc.files$merged$isoform), class.files$merged$isoform)
setdiff(class.files$merged$isoform, unique(junc.files$merged$isoform))

# output gffcomp_match
gffcomp_match_output <- gffcomp_match %>% 
  select(ref_gene_id, ref_id, qry_gene_id, qry_id) %>% 
  `colnames<-`(c("brain_associated_gene", "brain_isoform", "pancreas_associated_gene","pancreas_isoform")) %>% 
  mutate(brain_pancreas_common_isoform = paste0(brain_isoform,"_", pancreas_isoform))

## ---------- Output -----------------

write.table(junc.files$merged, paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted_junctions.txt"), quote = F, sep = "\t", row.names = F)
write.table(brain_pancreas_phenotype, paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/phenotype.txt"), quote = F, row.names = F, sep = "\t")
write.table(final_exp_fetal, paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted_expression.txt"), quote = F, row.names = F, sep = "\t")
write.table(merged_gtf,paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(class.files$merged,paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted_classification.txt"), sep = "\t", col.names = T, quote = F)
write.table(gffcomp_match_output, paste0(ROOT_DIR, "/9_tappAS_SK/1_Input/C_BrainvsPancreas/brain_pancreas_common_genes_isoforms.txt"), quote = F, row.names = F, col.names = T)




# check that all the isoforms are in the gff3 
gff3 <- read.table(paste0(ROOT_DIR,"10_characterisation_SK/Brain_vs_Pancreas/2_sqanti3/SFARI.gff3"), as.is = T, sep = "\t")
setdiff(rownames(final_exp),gff3$V1)
setdiff(gffcomp_match$qry_id,gff3$V1)
setdiff(gff3$V1,rownames(final_exp))
