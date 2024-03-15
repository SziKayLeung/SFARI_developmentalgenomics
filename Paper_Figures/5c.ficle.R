protein_coding_genes = read.table("/lustre/home/vc362/protein-coding-genes.txt")
FICLE <- list(

  # final transcript classifications of all merged samples
  Gene_class = lapply(list.files(path = dirnames$ficle, pattern = "final_transcript_classifications.csv", recursive = TRUE, full = T), 
                      function(x) {
                        print(x)
                        read.csv(x)}
                      )
  
)
names(FICLE$Gene_class) = lapply(list.files(path = dirnames$ficle, pattern = "final_transcript_classifications.csv", recursive = TRUE), 
                                    function(x) word(x, c(1), sep = fixed("/")))
Merged_gene_class_df <- all_summarise_gene_stats(Gene_class=FICLE$Gene_class,class_files=class.files$glob_SQ,Cpat$whole,Cpat$whole_noORF, names(FICLE$Gene_class))

#for(i in setdiff(protein_coding_genes,names(FICLE$Gene_class){
#  FICLE$Gene_class[[i]] = read.csv(paste0(i,"_final_transcript_classifications.csv"), recursive = TRUE, full = T))
#}

FICLE_class <- fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/TargetGenes/all_final_transcript_classifications.csv", data.table = F)
FICLE_class <- FICLE_class[!duplicated(FICLE_class), ]
FICLE_class <- FICLE_class %>% select(-isoform)
FICLE_class_final <- FICLE_class %>% filter(isoform %in% class.files$glob_SQ_annoGene$isoform)
protein_coding_genes_isoforms <- class.files$glob_SQ[class.files$glob_SQ$associated_gene %in% protein_coding_genes$V1,"isoform"]
setdiff(rownames(FICLE_class_final$isoform), protein_coding_genes_isoforms)
FICLENotProcessed <- setdiff(protein_coding_genes_isoforms, rownames(FICLE_class_final$isoform))
FICLEProcessedGenes <- unique(class.files$glob_SQ[class.files$glob_SQ$isoform %in% FICLE_class_final$isoform,"associated_gene"])
FICLENotProcessedGenes <- unique(class.files$glob_SQ[class.files$glob_SQ$isoform %in% FICLENotProcessed,"associated_gene"])
unique(setdiff(FICLENotProcessedGenes,FICLEProcessedGenes))

FICLE_class_final_exp <- merge(FICLE_class_final[,c("isoform","A5A3","ES","IR","NE_All")], class.files$glob_SQ[,c("isoform","preReads","postReads")], by = "isoform")
FICLE_class_final_exp <- merge(FICLE_class_final_exp, normWholeIsoform, by = "isoform", all.x = T)

FICLE_class_final_exp <- FICLE_class_final_exp %>% mutate(A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3),
                                                          ES = as.numeric(ES), ESpre = ES * normPre/sum(ES), ESpost = ES * normPost/sum(ES),
                                                          IR = as.numeric(IR), IRpre = IR * normPre/sum(IR), IRpost = IR * normPost/sum(IR),
                                                          A5A3 = as.numeric(A5A3), A5A3pre = A5A3 * normPre/sum(A5A3), A5A3post = A5A3 * normPost/sum(A5A3))

FICLEES <- FICLE_class_final_exp %>% select(isoform, ESpre, ESpost) %>% reshape2::melt(variable.name = "development", value.name = "normES") 
FICLEIR <- FICLE_class_final_exp %>% select(isoform, IRpre, IRpost) %>% reshape2::melt(variable.name = "development", value.name = "normIR") 
FICLEA5A3 <- FICLE_class_final_exp %>% select(isoform, A5A3pre, A5A3post) %>% reshape2::melt(variable.name = "development", value.name = "normA5A3") 

ggplot(FICLEES, aes(x = development, y = log10(normES))) + geom_boxplot()
ggplot(FICLEIR, aes(x = development, y = log10(normIR))) + geom_boxplot()
ggplot(FICLEA5A3, aes(x = development, y = log10(normA5A3))) + geom_boxplot()


t.test(normES ~ development, FICLEES)
t.test(normA5A3 ~ development, FICLEA5A3)
t.test(normIR ~ development, FICLEIR)

row.names(FICLE_class) <- FICLE_class$isoform
FICLE_class <- FICLE_class %>% mutate_if(is.character, as.numeric)
FICLE_class <- FICLE_class[complete.cases(FICLE_class), ]


count = 1
AS_events <- data.frame()
for(i in colnames(FICLE_class)){
  AS_events[count,1] <- i
  AS_events[count,2] <- sum(FICLE_class[[i]])
  count <- count + 1
}
colnames(AS_events) <- c("AS_event","Number")
ggplot(AS_events, aes(x = reorder(AS_event,-Number), y = Number)) + geom_bar(stat = "identity")


plot_dendro_Tgene(dirnames$ficle,"DLGAP5")
DLGAP5ES <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"ES_events_counts.csv") %>% arrange(-numEvents) 
DLGAP5NE <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"NE_transcript_counts.csv") %>% arrange(-numNovelExons) 
DLGAP5IR <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"IR_events_counts.csv") %>% arrange(-numEvents) 

table(FICLE_class_final_exp)


DLGAP5Prenatal <- class.files$glob_SQ[class.files$glob_SQ$associated_gene == "DLGAP5" & class.files$glob_SQ$DevStatus == "prenatal",] %>% arrange(-preReads)
DLGAP5Iso <- data.frame(
  Isoform = unlist(DLGAP5Iso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "DLGAP5","transcript_id"])),
    `Pre-Top` = DLGAP5Prenatal[["isoform"]][1:10],
    `Pre-ES` = as.character(DLGAP5ES[,c("transcriptID")][1:5]),
    `Pre-NE` = as.character(DLGAP5NE[,c("transcriptID")][1:5]),
    `Pre-IR` = as.character(DLGAP5IR$transcriptID)
  )),
  Category = rep(names(DLGAP5Iso), lengths(DLGAP5Iso))
)
DLGAP5Iso$colour <- NA

ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
            isoList = c(as.character(DLGAP5Iso$Isoform)),
            selfDf = DLGAP5Iso, gene = "DLGAP5")

PKM <- list(
  IF = plotIFWholebyGene("PKM",dirnames$DIU),
  GeneExp = plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$whole_group,gene="PKM"),
  T1Exp = plot_trans_exp_lifetime("ONT15_1709_8230",class.files$glob_targ_SQ,Exp$whole_group),
  T2Exp = plot_trans_exp_lifetime("ONT15_1709_8237",class.files$glob_targ_SQ,Exp$whole_group),
  Track = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
              isoList = c("ONT15_1709_8230","ONT15_1709_8237",RefIsoforms$PKM$transcript_id),
              colours = c(wes_palette("Cavalcanti1")[5],wes_palette("Cavalcanti1")[2],"#0C0C78"),
              simple=TRUE)
)

pdf("PKM_DIU.pdf", width = 12, height = 14)
plot_grid(
  plot_grid(PKM$IF,PKM$GeneExp, labels = c("i","ii")),
  plot_grid(PKM$Track, labels = c("iii")),
  plot_grid(PKM$T1Exp,PKM$T2Exp, labels = c("iv","v")),
  ncol = 1, rel_heights = c(0.5,0.1,0.3)
)
dev.off()
