#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: R.Bamford, K.Chundru, S.Leung (S.K.Leung@exeter.ac.uk)
## output a list of transcripts that are documented in the TALON gtf but not in the abundance file
## sanity check for TALON pipeline
## Input:
##  --gtf = TALON gtf generated from talon_create_GTF
##  --counts = TALON abundance file
##  --output = output path and name of file
## --------------------------------

## ---------- packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(optparse))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-p", "--phenotype"), type="character", default=NULL, 
              help="Phenotype file", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="Expression file"),
  make_option(c("-c", "--classfile"), type="character", default=NULL, 
              help="classification file"),
  make_option(c("-t", "--threshold"), type="numeric", default=10L, 
              help="read threshold for DESeq"),
  make_option(c("-d", "--design"), type="character", default=NULL, 
              help="<development><sex>"),
  make_option(c("-l", "--level"), type="character", default="transcript", 
              help="<transcript><gene>"),
  make_option(c("-i", "--dataset"), type="character", default="whole", 
              help="<targeted><whole>"),
  make_option(c("-o", "--directory"), type="character", 
              help="output directory")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ---------- Input files -----------------

# phenotype
message("Read in: ", opt$phenotype)
phenotype <- read.csv(opt$phenotype) 
nrow(phenotype)
# remove early postnatal samples
removePhenotype <- phenotype[phenotype$group == "Postnatal" & phenotype$age < 1,]
message("Removing samples", removePhenotype$sample)
phenotype <- phenotype[!phenotype$sample %in% removePhenotype$sample,]
nrow(phenotype)

# expression (demux csv - make sure first column is isoform; all samples are there and that the sample names match up with those in the phenotype file)
message("Read in: ", opt$expression)
expression <- fread(opt$expression, data.table=F)

# classification file
class.files <- fread(opt$classfile, sep="\t", header=T, data.table=F)
row.names(class.files) <- class.files$isoform


## ---------- data wrangle -----------------

## --- class.files --- 
# column (filter_result) present in classification file, make sure to only select true isoforms 
# otherwise assumed input classification file has been filtered prior 
if(exists("filter_result", where = class.files) == TRUE){
  class.files <- class.files[which(class.files$filter_result=="Isoform"),]
}
class.files <- class.files[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory","exons")]

## --- expression --- 
# subset samples based on the phenotype file
expression <- expression[,c(1,which(colnames(expression)%in%phenotype$sample))]

# subset expression to only the transcripts in classification file
expression <-  right_join(expression,class.files[,c("isoform","associated_gene")], by = "isoform")

# filter rows with 0 reads
expression <- expression[which(rowSums(expression[,2:(ncol(expression)-1)])>0),]

# remove associated gene column
if(opt$level == "gene"){
  expression <- as.data.frame(expression %>%
                                group_by(associated_gene) %>%
                                select(-isoform) %>% summarize_all(.funs = list(sum)))
  rownames(expression)<- expression$associated_gene
  expression <- expression %>% dplyr::select(-associated_gene)
  
}else{
  rownames(expression) <- expression$isoform
  expression <- expression %>% dplyr::select(-associated_gene, isoform)
}

## --- phenotype --- 
# factor group level 
if(opt$design == "development"){
  phenotype$group <- as.factor(phenotype$group)
  phenotype$group <- relevel(phenotype$group,"Postnatal")
  str(phenotype$group)
}else{
  phenotype$sex <- as.factor(phenotype$sex)
  phenotype$sex <- relevel(phenotype$sex,"M")
  str(phenotype$sex)
}


## matching samples
if(nrow(phenotype) != ncol(expression)){
  message("Mismatched number of samples in phenotype and expression file")
  message("Using common samples in both files")
  col_match <- intersect(phenotype$sample, colnames(expression))
  cat("Number of samples:", length(col_match),"\n")
  phenotype <- phenotype %>% filter(sample %in% col_match) 
  expression <- expression %>% dplyr::select(phenotype$sample)
}


## ---------- run DESeq -----------------

if(opt$design == "development"){
  message("Design: ~ group")
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), colData = phenotype, design = ~ group)
}else{
  message("Design: ~ sex")
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), colData = phenotype, design = ~ sex)
}
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))

cat("Number of isoforms before filtering:", nrow(dds),"\n")
cat("Filtering isoform on count threshold:", opt$threshold,"\n")
dds <- dds[rowSums(counts(dds)) > opt$threshold, ] 
cat("Number of isoforms after filtering:", nrow(dds),"\n")

cat("Running Wald test\n")
dds_Wald <- DESeq(dds, test="Wald")
if(opt$level == "gene"){
  res_Wald <- as.data.frame(results(dds_Wald)) %>% tibble::rownames_to_column("associated_gene") %>% arrange(padj) 
}else{
  res_Wald <- as.data.frame(results(dds_Wald)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj) 
}

deseq_output <- list(dds_Wald, res_Wald, norm)
names(deseq_output) <- c("dds_Wald", "res_Wald","norm_counts")

## ---------- output -----------------

# wald results and annotate with classification file
if(opt$level == "gene"){
  deseq_output$norm <- deseq_output$norm %>% dplyr::rename("associated_gene" = "isoform")
  deseq_output$anno_res <- deseq_output$res_Wald
  deseq_output$anno_norm <- deseq_output$norm
}else{
  deseq_output$anno_res <-  right_join(deseq_output$res_Wald,class.files, by = "isoform")
  deseq_output$anno_norm <- right_join(deseq_output$norm,class.files, by = "isoform")
}

fwrite(deseq_output$anno_res, file=paste0(opt$directory, "/DESeq2_", opt$dataset,"_",opt$level,"_",opt$design,"_resAll.csv"))
fwrite(deseq_output$anno_res[which(deseq_output$anno_res$padj<0.05),], file=paste0(opt$directory, "/DESeq2_", opt$dataset,"_", opt$level, "_", opt$design,"_resSig.csv"))

# normalisation
fwrite(deseq_output$anno_norm,file=paste0(opt$directory, "/DESeq2_", opt$dataset,"_", opt$design,"_normAll.csv"))
normSig <- deseq_output$anno_norm[which(deseq_output$anno_norm$isoform%in%deseq_output$res_Wald[which(deseq_output$res_Wald$padj<0.05),]$isoform),]
fwrite(normSig,file=paste0(opt$directory, "/DESeq2_", opt$dataset,"_", opt$design,"_normSig.csv"))
