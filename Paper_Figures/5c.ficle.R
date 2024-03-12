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

plot_dendro_Tgene(dirnames$ficle,"DLGAP5")
DLGAP5ES <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"ES_events_counts.csv") %>% arrange(-numEvents) 
DLGAP5NE <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"NE_transcript_counts.csv") %>% arrange(-numNovelExons) 
DLGAP5IR <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"IR_events_counts.csv") %>% arrange(-numEvents) 


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
