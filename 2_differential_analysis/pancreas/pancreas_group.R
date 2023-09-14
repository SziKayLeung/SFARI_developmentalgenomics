## ---------- Script -----------------
##
## Purpose: perform differential analysis on pancreas targeted datasets using linear regression
## Transcript level separate analysis
## Gene level separate analysis
##
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
##
## --------- Notes -----------------
## Expression


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

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
run_DESeq2 <- function(expression,phenotype,exprowname,threshold=10,controlname="Control",design="time_series",interaction="Off",test){
  
  
  # input phenoptype characteristics as factor
  phenotype$group <- as.numeric(phenotype$group)
  if(design == "time_series"){
    phenotype$time <- as.factor(phenotype$time)
  }
  
  # common columns between phenotype and expression
  col_match <- intersect(phenotype$sample, colnames(expression))
  cat("Number of samples:", length(col_match),"\n")
  
  # ensure expression column and phenotype rows are in the same order
  rownames(expression) <- expression[[exprowname]]
  phenotype <- phenotype %>% filter(sample %in% col_match) 
  expression <- expression %>% dplyr::select(phenotype$sample)
  rownames(phenotype) <- phenotype$sample
  phenotype <- phenotype %>% select(-sample)
  if(all(colnames(expression) == rownames(phenotype))==FALSE){
    print("ERROR: rownames and colnames in expression matrix and phenotype are not in the same order")
  }
  
  if(nrow(phenotype) != ncol(expression)){
    print("Mismatched number of samples in phenotype and expression file")
  }
  
  cat("Expression file:\n")
  print(head(expression))
  # create DESeq2 data object
  # normalises reads within DESeq2
  # counts need to be rounded to the nearest integer
  # design is linear model – last expression is tested for differential expression
  
  
  if(design == "time_series"){
    cat("Design: time_series\n")
    if(interaction == "Off"){
      cat("Modelling interaction effect: off\n")
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                    colData = phenotype,
                                    design = ~ group + time)
      
    }else{
      cat("Modelling interaction effect: on\n")
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                    colData = phenotype,
                                    design = ~ group + time + group:time)
      
    }
  }else{
    cat("Design: case_control\n")
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                  colData = phenotype,
                                  design = ~ group)
  }
  
  
  
  # estimate size factors to account for differences in sequencing depth
  # if all samples have exact same sequencing depth, size factor should be ~ 1
  dds <- estimateSizeFactors(dds)
  
  # pre-filtering the data set
  # remove the rows of the DESeqDataSet that have no counts, or only a single count across all samples
  # arbitrary threshold = 10 (default)
  # minimum 2 samples per group, and expect minimum 2 FL reads per isoform
  cat("No of isoforms before filtering:", nrow(dds),"\n")
  cat("Filtering isoform on count threshold (minimum):", threshold,"\n")
  dds <- dds[rowSums(counts(dds)) >= threshold, ]
  cat("No of isoforms after filtering:", nrow(dds),"\n")
  
  #p <- cbind(reshape2::melt(sizeFactors(dds)), reshape2::melt(colSums(counts(dds)))) %>%
  #  `colnames<-`(c("sizefactors", "nreads")) %>%
  #  tibble::rownames_to_column("sample") %>%
  #  mutate(sample = str_remove(sample,"ont_")) %>%
  #  ggplot(., aes(x = sizefactors, y = nreads)) + geom_point() +
  #  geom_label_repel(aes(label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  #  theme_bw() + labs(y = "Number of reads", x = "Size Factors")
  
  
  # normalization to stabilize variance (regularized logarithm)
  #rld <- rlog(dds, blind = FALSE)
  
  # PCA plot
  #pcaData <- plotPCA(rld, intgroup = c("time", "group"), returnData = TRUE)
  #pcaData$sample <- sapply(pcaData$name, function(x) str_remove(x, "ont_"))
  #percentVar <- round(100 * attr(pcaData, "percentVar"))
  #p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape = group.1)) +
  #  geom_point(size =3) +
  #  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance")) +
  #  coord_fixed() +
  #  geom_label_repel(aes(label = sample), segment.color = 'grey50')
  
  
  # normalised counts for downstream plotting
  norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))
  
  # run differential analysis
  if(test == "Wald"){
    cat("Running Wald test\n")
    #colData(dds)$group <- relevel(colData(dds)$group, controlname)
    dds_output <- DESeq(dds, test="Wald")
    output_names <- c("dds_Wald", "res_Wald","norm_counts", "stats_Wald")
  }else if(test == "LRT"){
    cat("Running LRT test\n")
    # reduced model removes the interaction term 
    # significant DE genes will represent those genes that have differences in the effect of genotype over time
    dds_output <- DESeq(dds, reduced=~group+time, test="LRT")
    output_names <- c("dds_LRT", "res_LRT","norm_counts", "stats_LRT")
  }else{
    print("test either <Wald/LRT>")
  }
  
  res <- as.data.frame(results(dds_output)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj)
  stats <- as.data.frame(mcols(dds_output))
  output <- list(dds_output, res, norm, stats)
  names(output) <- output_names
  
  return(output)
}

# plots
# only keep the isoforms that are in the sqanti classification file and not in the expression file
anno_DESeq2 <- function(deseq_output,class.files,phenotype,controlname="Control",level,sig=0.05){
  
  #phenotype$group <- as.factor(phenotype$group)
  #phenotype$group <- relevel(phenotype$group,controlname)
  
  if(level == "transcript"){
    print("Processing via transcript")
    deseq_output$anno_res <-  merge(deseq_output[[2]],class.files[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory")], by = "isoform", all.x = T)
    # normalised counts
    deseq_output[[3]] <- merge(deseq_output[[3]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
    deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)
  }else{
    print("Processing via gene")
    class.files <- class.files %>% mutate(PB_associated_gene = word(isoform,c(2), sep = fixed(".")))
    gene_class_files <- unique(class.files[,c("associated_gene","PB_associated_gene")])
    deseq_output$anno_res <- merge(deseq_output[[2]], gene_class_files, by.x = "isoform", by.y = "PB_associated_gene")
    
    # normalised counts
    deseq_output[[3]] <- merge(deseq_output[[3]],gene_class_files, by.x = "isoform", by.y = "PB_associated_gene", all.x = T)
    deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)
    
    deseq_output[[4]] <- merge(deseq_output[[4]], gene_class_files, by.x = 0, by.y = "PB_associated_gene")
  }
  
  deseq_output$anno_res <- deseq_output$anno_res %>% filter(padj < sig) %>% arrange(padj)
  
  return(deseq_output)
}



## ---------- input -----------------

# directory names
dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/",
  combined =  "/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Figures/"
)

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$root, "ONT_pheno.txt"), 
  expression = paste0(dirnames$combined, "demux_fl_count.csv"),
  classfiles = paste0(dirnames$combined, "SQANTI3_whole_RulesFilter_result_classification.targetgenes_counts_filtered.txt")
)


input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = "\t", header = T) %>% select(sample,time,group)
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$expression <- read.csv(input_files$expression) %>% .[.$id != "0",] %>% filter(id %in% input$classfiles$isoform)

input$ontPhenotype <- input$ontPhenotype %>% select(sample, time)
colnames(input$ontPhenotype) <- c("sample","group")

## ---------- ONT: Creating DESeq2 object and analysis -----------------

# run DESeq2
ontResTran <- list(
  wald = run_DESeq2(test="Wald",input$expression,input$ontPhenotype,threshold=10,exprowname="id",controlname="M",design="case_control",interaction="On")
)

ontResTranAnno <- lapply(ontResTran, function(x) anno_DESeq2(x,input$classfiles,input$ontPhenotype,controlname="M",level="transcript",sig=0.1))


## ---------- Output -----------------

saveRDS(ontResTranAnno, file = paste0(dirnames$output, "/Ont_DESeq2TranscriptLevel.RDS"))