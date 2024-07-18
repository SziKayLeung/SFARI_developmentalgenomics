##Szi Kay run_DESSeq2.R script


library(dplyr) #needed for %>% command
library(DESeq2)
library(ggrepel)
library(data.table)
source("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/run_DESeq2.R")

#Phenotype file
phenotype=read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/00_metadata/WholeTargetedphenotype_fixedsex.csv")
phenotype=phenotype[which(startsWith(phenotype$sample, 'Tar')),]

#Expression file (demux csv - make sure first column is isoform; all samples are there and that the sample names match up with those in the phenotype file)
expression=fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/WholeTargeted_demux.csv",data.table=F)

#selecting only targeted genes
#expression = expression[,1:41]
expression = expression[,c(1,which(colnames(expression)%in%phenotype$sample))]

#adding 1 to each cell
#expression <- expression %>% mutate(across(where(is.numeric), ~ .x +1))

#Filter expression file for targeted genes (but not filtered on abundance)
#input.class.files="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/1_SQANTI3/SQANTI3_collapse_options_RulesFilter_result_classification.targetgenes_counts.txt"

input.class.files="/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/unfiltered/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification.txt"

class.files=fread(input.class.files,sep="\t",header=T,data.table=F)
#print(head(class.files))
class.files <- class.files[which(class.files$filter_result=="Isoform"),]

expression <-  right_join(expression,class.files[,c("isoform","associated_gene")], by = "isoform")# %>% dplyr::select(-associated_gene)
#expression = expression[which(rowSums(expression[,2:41])>0),]
expression <- expression[which(rowSums(expression[,2:(ncol(expression)-1)])>0),]
tgenes<-read.table('/gpfs/mrc0/projects/Research_Project-MRC148213/vc362/fetalBrain/genes.txt', stringsAsFactors=F)
expression <- expression[which(expression$associated_gene%in%tgenes$V1),]
expression <- expression %>% dplyr::select(-associated_gene)

#fwrite(expression,file="wholetargetedfiltered_demux.csv", quote=F)

#nThread=16 helps with big files
#print(head(expression))
phenotype$sex <- as.factor(phenotype$sex)
#phenotype$group <- as.factor(phenotype$group)
#phenotype$group <- factor(phenotype$group, levels=c("Early", "Mid...)

#phenotype$group <- relevel(phenotype$group,"NeuN")
phenotype$sex <- relevel(phenotype$sex,"M")

#phenotype$group <- relevel(phenotype$group,"Postnatal")
#phenotype$sex <- relevel(phenotype$sex,"M")

#str(phenotype$group)
str(phenotype$sex)

col_match <- intersect(phenotype$sample, colnames(expression))
cat("Number of samples:", length(col_match),"\n")


#expression <- expression %>% tibble::column_to_rownames("isoform")
rownames(expression)<-expression$isoform
phenotype <- phenotype %>% filter(sample %in% col_match) 
expression <- expression %>% dplyr::select(phenotype$sample)

if(nrow(phenotype) != ncol(expression)){
	print("Mismatched number of samples in phenotype and expression file")
}

cat("Expression file:\n")


dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), colData = phenotype, design = ~ sex)

dds <- estimateSizeFactors(dds)

norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))
phenotype=phenotype[which((phenotype$group=='Prenatal'&phenotype$age>30)==F & (phenotype$group=='Postnatal'&phenotype$age<2)==F),]

input.class.files="/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt"
class.files=read.table(input.class.files,sep="\t",as.is=T,header=T)

#dds<-dds[which(rownames(dds)%in%class.files$isoform),which(dds$sample%in%phenotype$sample)]
dds<-dds[,which(dds$sample%in%phenotype$sample)]

threshold=10

cat("Number of isoforms before filtering:", nrow(dds),"\n")
cat("Filtering isoform on count threshold:", threshold,"\n")
dds <- dds[rowSums(counts(dds)) > threshold, ]
#dds <- dds[which(rowSums(counts(dds)>0) > 1), ]
cat("Number of isoforms after filtering:", nrow(dds),"\n")

cat("Running Wald test\n")
dds_Wald <- DESeq(dds, test="Wald")
res_Wald <- as.data.frame(results(dds_Wald)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj) 

deseq_output <- list(dds_Wald, res_Wald, norm)

names(deseq_output) <- c("dds_Wald", "res_Wald","norm_counts")

fwrite(deseq_output$res_Wald, file="output_targeted_sex_removinglateprenatal.csv")
fwrite(deseq_output$res_Wald[which(deseq_output$res_Wald$padj<0.05),], file="filtered_output_targeted_sex_removinglateprenatal.csv")

fwrite(deseq_output$norm,file="output_norm_targeted_sex_removinglateprenatal.csv")
fwrite(deseq_output$norm[which(deseq_output$norm$isoform%in%deseq_output$res_Wald[which(deseq_output$res_Wald$padj<0.05),]$isoform),],file="filtered_output_norm_targeted_sex_removinglateprenatal.csv")

##Targeted

deseq_output$anno_res <-  right_join(deseq_output$res_Wald,class.files[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory")], by = "isoform")

fwrite(deseq_output$anno_res,file="anno_targeted_sex_removinglateprenatal.csv", quote=F)

##deseq_output[[3]] <- merge(deseq_output[[3]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
##deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)

##write.csv(deseq_output[[3]],file="anno_output.csv")

SFARI=fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/12_deseq2/SFARI-Gene_genes_01-23-2023release_03-24-2023export.csv",data.table=F)

##class.files=read.table(SFARI,sep="\t",as.is=T,header=T)

SFARI_hits <-  right_join(deseq_output$anno_res,SFARI[,c("gene-symbol","associated_gene")], by = "associated_gene")

fwrite(SFARI_hits,file="Targeted_SFARI_hits_sex_removinglateprenatal.csv", quote=F)
