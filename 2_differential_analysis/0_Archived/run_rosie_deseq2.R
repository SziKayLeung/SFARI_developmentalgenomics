##Szi Kay run_DESSeq2.R script


suppressMessages(library(dplyr)) 
suppressMessages(library(DESeq2))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
output_dir="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/9_tappAS_SK"
source("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/run_DESeq2.R")
phenotype=read.csv("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/phenotype.csv")
expression=fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/demux_fl_count.csv",nThread=16)

print(head(expression))
phenotype$time <- as.factor(phenotype$time)
phenotype$group <- as.factor(phenotype$group)
phenotype$group <- relevel(phenotype$group,"M")
phenotype$time <- relevel(phenotype$time,"Prenatal")

str(phenotype$group)
str(phenotype$time)

col_match <- intersect(phenotype$sample, colnames(expression))
cat("Number of samples:", length(col_match),"\n")

expression <- expression %>% tibble::column_to_rownames("isoform")
#rownames(expression) <- expression[["isoform"]]
phenotype <- phenotype %>% filter(sample %in% col_match) 
expression <- expression %>% dplyr::select(phenotype$sample)

if(nrow(phenotype) != ncol(expression)){
  print("Mismatched number of samples in phenotype and expression file")
}

#cat("Expression file:\n")


dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), colData = phenotype, design = ~ group + time)

dds <- estimateSizeFactors(dds)

threshold=10

cat("No of isoforms before filtering:", nrow(dds),"\n")
cat("Filtering isoform on count threshold:", threshold,"\n")
dds <- dds[rowSums(counts(dds)) > threshold, ] 
cat("No of isoforms after filtering:", nrow(dds),"\n")

cat("Running Wald test\n")
dds_Wald <- DESeq(dds, test="Wald")
res_Wald <- as.data.frame(results(dds_Wald)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj) 
norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))

output <- list(dds_Wald, res_Wald, norm)
names(output) <- c("dds_Wald", "res_Wald","norm_counts")
write.csv(output$res_Wald, file=paste0(output_dir,"/deseq2_res_wald.csv"))
write.csv(output$norm, file=paste0(output_dir,"/deseq2_norm.csv"))

input.class.files="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/11_cupcake/merged_collapse_classification.txt"
class.files=read.table(input.class.files,sep="\t",as.is=T,header=T)

# top 50 hits
output_50 <- output
output_50$res_Wald <- output_50$res_Wald[1:50,]
output_50$norm <- output_50$norm %>% filter(isoform %in% output_50$res_Wald$isoform)
output_50$class.files <- class.files %>% filter(isoform %in% output_50$res_Wald$isoform)
annoOutput <- list()

write.table(output_50$res_Wald$isoform, file=paste0(output_dir,"/deseq2_50sig_isoform.tsv"),quote=F,row.names = F,col.names = F)

merge(output_50$norm, phenotype)

annoOutput = anno_DESeq2(output_50,class.files,phenotype,controlname="F",level="transcript",sig=0.05)

