# global parameters
# need disease_list.csv for AD genes, SZ genes 

# updated autism SFARI list 
diseasegenelists <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/0_utilities/disease_list/"
SFARI <- read.csv(paste0(diseasegenelists,"SFARI-Gene_genes_07-17-2023release_09-26-2023export.csv"),header=T)
SFARI_CLASS_1_2_S <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2 |SFARI$syndromic == 1)
SFARI_CLASS_1_2 <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2)
disease_list <- list(
  SCHEMA = read.table(paste0(diseasegenelists,"SCHEMA_Oct2023.csv"), sep=",", header=T, stringsAsFactors = F),
  DDG2P = read.table(paste0(diseasegenelists,"DDG2P_26_9_2023.csv"), sep=",", header=T, stringsAsFactors = F) %>% filter(confidence.category != "limited")
)


# 1. Datawrangle the disease list  
# aim: remove white space or empty records from each column of the DiseaseList.csv and save genes as new list 
# output: rehash_disease_list  
rehash_disease_list <- list(
  SZ = unique(as.factor(disease_list$SCHEMA$Gene)),
  Autism = SFARI_CLASS_1_2,
  DDG2P = unique(as.factor(disease_list$DDG2P$gene.symbol))
)

# 2. Number of transcripts and genes in disease
disease_all <- function(input.class.files){
  disease <- unique(unlist(rehash_disease_list))
  alldisease_class.files <- input.class.files[input.class.files$associated_gene %in% disease,]
  
  totaltranscripts <- nrow(alldisease_class.files)
  totalgenes <- length(unique(alldisease_class.files$associated_gene))
  cat("Number of transcripts to disease genes:",totaltranscripts,"\n")
  cat("Number of genes to disease genes:",totalgenes,"\n")
  
  # more than one isoform
  morethan1 <- alldisease_class.files %>% group_by(associated_gene) %>% tally() %>% filter(n > 1) %>% nrow()
  cat("Number of disease genes with more than isoform:",morethan1,"(",round(morethan1/totalgenes*100,2),"%)","\n")
  
  # novel transcripts 
  novel <- alldisease_class.files %>% filter(associated_transcript == "novel") %>% nrow()
  cat("Number of novel transcripts of disease genes:",novel,"(",round(novel/totaltranscripts*100,2),"%)","\n")
}



# 3. Subset classification file of just the genes interested for each disease 
# aim: Loop through the rehash_disease_list and extract from classification file the rows mapped to those genes 
# input: rehash_disease_list 
# output: disease_class.files 
subset_disease <- function(input.class.files){
  # use count to loop and save each diseaes gene list into new list entry 
  disease_class.files <- list()
  count=1
  for(i in rehash_disease_list){
    disease_class.files[[count]] <- input.class.files[input.class.files$associated_gene %in% i,]
    disease_class.files[[count]]$disease <- i[1]
    count = count + 1
  }
  names(disease_class.files) <- names(rehash_disease_list)
  return(disease_class.files)
}



# 4. General Basic Stats about genes associated with disease 
# aim: Extract the number of associated isoforms, novel and known and the differen SQANTI2 categories 
# input: disease_class.files 
# output: output basic stats
tabulating_sqanti_num <- function(type_class_file){
  dat <- data.frame()
  count=1
  for(i in type_class_file){
    # Total Unique Isoforms: tabulated by number of rows
    isoforms <- dim(i)[1]
    # Total Unique Genes: Remove novel genes and count 
    annotated_genes <- nrow(unique(i[!grepl("NOVELGENE",i$associated_gene),"associated_gene"]))
    # Number of Annotated Isoforms (FSM, ISM) 
    annotated_isoforms <- paste0(nrow(i[i$associated_transcript != "novel",])," (",
                                 round(nrow(i[i$associated_transcript != "novel",])/isoforms * 100,2),"%)")
    novel_isoforms <- paste0(nrow(i[i$associated_transcript == "novel",])," (",
                             round(nrow(i[i$associated_transcript == "novel",])/isoforms * 100,2),"%)")
    
    # 9 levels of structural cateogory 
    struct <- vector("numeric", 9)
    for(num in 1:length(levels(i$structural_category))){
      struct[num] <- nrow(i[i$structural_category == levels(i$structural_category)[num],]) 
    }
    
    dat[1:4,count] <- rbind(annotated_genes, isoforms, annotated_isoforms, novel_isoforms)
    dat[5:13,count] <- struct
    colnames(dat)[count] <- names(type_class_file)[count]
    count = count + 1
  }
  row.names(dat) <- append(c("Total Number of Detected Genes", "Total Number of Isoforms", 
                             "Number and % of Annotated Isoforms", "Number and % of Novel Isoforms"),
                           levels(i$structural_category))
  
  dat
}

num_DTE <- function(deseq2_list, gene, colname){
  print(gene)
  output = data.frame(
    gene = gene, 
    number = nrow(subset(deseq2_list, associated_gene == "gene")) 
  )
  output = output %>% dplyr::rename(!!colname := number)
  return(output)
}

# 5. IR vs NMD vs NMD_IR
# aim: Tabulating for each gene in the disease list the number of IR, NMD and Fusion transcripts 
# input: rehash_disease_list 
# output: number of associated IR vs NMD transcripts for each gene: disease_num
num_disease <- function(disease_list, input.class.files){
  output <- data.frame()
  for(i in 1:length(disease_list)){
    df <- input.class.files[input.class.files$associated_gene == disease_list[i],]
    output[i,1] <- disease_list[i]
    output[i,2] <- nrow(df)
    # Number of Annotated Isoforms (FSM, ISM) 
    output[i,3] <- paste0(nrow(df[df$associated_transcript != "novel",])," (",
                                 round(nrow(df[df$associated_transcript != "novel",])/nrow(df) * 100,2),"%)")
    output[i,4] <- paste0(nrow(df[df$associated_transcript == "novel",])," (",
                          round(nrow(df[df$associated_transcript == "novel",])/nrow(df) * 100,2),"%)")
    
    # isoform, exon length
    output[i,5] <- paste0("median=", median(df$length)," (SD=",
                           round(sd(df$length),0),"), range=", min(df$length),"-",max(df$length))
    output[i,6] <- paste0("median=", median(df$exons)," (SD=",
                           round(sd(df$exons),0),"), range=", min(df$exons),"-",max(df$exons))
    
    # 9 levels of structural cateogory 
    struct <- vector("numeric", 9)
    for(num in 1:length(levels(df$structural_category))){
      struct[num] <- nrow(df[df$structural_category == levels(df$structural_category)[num],]) 
    }
    
    output[i,7:15] <- struct
    
  }
  
  colnames(output) <- append(c("associated_gene", "Total Number of Isoforms", 
                             "Number and % of Annotated Isoforms", "Number and % of Novel Isoforms",
                             "Isoform lengths","Isoform Exons"),
                           levels(df$structural_category))
  full_list <- list(output)
  disease_num <- Reduce(function(...) merge(..., by='associated_gene', all.x=TRUE), full_list) 
  return(disease_num)
}


### Apply all functions into one 

# 1. All disease
disease_all(input.class.files)

# 2. Subset classification file of just the genes interested for each disease 
disease_class.files <- subset_disease(input.class.files)

# 3. General Basic Stats about genes associated with disease 
disease_stats <- tabulating_sqanti_num(disease_class.files)

# 4. Tabulate info at gene level associated with disease
num_disease_output_all <- lapply(rehash_disease_list, function(x) num_disease(x, input.class.files))

#write.csv(do.call(rbind, num_disease_output_all) %>% mutate(disease = word(rownames(.),c(1),  sep = fixed ('.'))) %>% distinct(),
#          paste0(dirnames$output, "DiseaseStats.csv"))

diseaseTargeted <- list()

diseaseTargeted$SZSex <- do.call(rbind, lapply(rehash_disease_list$SZ, function(x) 
  num_DTE(TargetedDESeq2Sig$sex, x, "TargetedDTESex")))

diseaseTargeted$SZGroup <- do.call(rbind, lapply(rehash_disease_list$SZ, 
                                                  function(x) num_DTE(TargetedDESeq2Sig$age, x, "TargetedDTEGroup")))

diseaseTargeted$AutismSex <- do.call(rbind, lapply(rehash_disease_list$Autism, function(x) 
  num_DTE(TargetedDESeq2Sig$sex, x, "TargetedDTESex")))

diseaseTargeted$AutismGroup <- do.call(rbind, lapply(rehash_disease_list$Autism, 
                                                  function(x) num_DTE(TargetedDESeq2Sig$age, x, "TargetedDTEGroup")))

saveRDS(diseaseTargeted, file = "diseaseTargeted.rds")

diseaseTargeted$DDG2PSex <- do.call(rbind, lapply(rehash_disease_list$DDG2P, function(x) 
  num_DTE(TargetedDESeq2Sig$sex, x, "TargetedDTESex")))

diseaseTargeted$DDG2PGroup <- do.call(rbind, lapply(rehash_disease_list$DDG2P, 
                                                     function(x) num_DTE(TargetedDESeq2Sig$age, x, "TargetedDTEGroup")))

