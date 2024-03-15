
message("Number of RNA transcripts protein coding: ", length(unique(proteinInput$t2p.collapse$pb_accs)))
message("Number of RNA isoforms: ", length(unique(proteinInput$t2p.collapse$base_acc)))

SC_ROOT <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation"
load(file = paste0(SC_ROOT,"/proteinInputWhole.RData"))
proteinCollapse <- readRDS("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/IsoProtein_Whole_DESeq2TranscriptLevel.RDS")
proteinCollapse$pWald$anno_res

diffTranscripts <- unique(proteinInput$t.class.files[proteinInput$t.class.files$isoform %in% WholeDESeqSig$age$isoform,"base_acc"])

nrow(WholeDESeqSig$age)
nrow(proteinCollapse$pWald$anno_res)
nrow(proteinCollapse$pWald$anno_res %>% filter(isoform %in% diffTranscripts))
View(proteinCollapse$pWald$anno_res %>% filter(isoform %in% diffTranscripts))

root <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/6a_longReadProteogenomics/5_calledOrfs"
input$peptide_orf <- as.data.frame(rtracklayer::import("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/6a_longReadProteogenomics/5_calledOrfs/Whole.gtf"))
input$peptide_orf$isoform <- gsub("*\\_ORF_[0-9]*","",input$peptide_orf$transcript_id)

input$merged_peptide_gtf <- rbind(gtf$glob_targ[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                                  input$peptide_orf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])



suppressMessages(library("cowplot"))
redundant_iso = proteinInput$t.class.files[proteinInput$t.class.files$base_acc  == "ONT20_1826_53630","isoform"]
redundant_MBP_iso = proteinInput$t.class.files[proteinInput$t.class.files$base_acc  == "ONT18_5132_10030","isoform"]
redundant_CAMK2A_iso = proteinInput$t.class.files[proteinInput$t.class.files$base_acc  == "ONT5_8947_2464","isoform"]



NNATIsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    `RNA Transcript` = redundant_iso[1:30],
    `RNA Isoform` = as.character(unique(input$peptide_orf[input$peptide_orf$isoform %in% redundant_iso[1:30],"transcript_id"]))
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
NNATIsoList$colour <- c(rep(NA,length(NNATIsoList$Category[NNATIsoList$Category != "DTE"])))

CAMK2AIsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    `RNA Transcript` = redundant_CAMK2A_iso[1:30],
    `RNA Isoform` = as.character(unique(input$peptide_orf[input$peptide_orf$isoform %in% redundant_CAMK2A_iso[1:30],"transcript_id"]))
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
CAMK2AIsoList$colour <- c(rep(NA,length(CAMK2AIsoList$Category[CAMK2AIsoList$Category != "DTE"])))
ggTranPlots(inputgtf=input$merged_peptide_gtf, classfiles=proteinInput$t.class.files, isoList = c(as.character(CAMK2AIsoList$Isoform)), selfDf = CAMK2AIsoList)

MBPIsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    `RNA Transcript` = c("ONT18_5132_2313", setdiff(redundant_MBP_iso[1:50],"ONT18_5132_2313")),
    `RNA Isoform` = as.character(unique(input$peptide_orf[input$peptide_orf$isoform %in% c("ONT18_5132_2313", "ONT18_5132_10030"),"transcript_id"]))
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
MBPIsoList$colour <- c(rep(NA,length(MBPIsoList$Category[MBPIsoList$Category != "DTE"])))
ggTranPlots(inputgtf=input$merged_peptide_gtf, classfiles=proteinInput$t.class.files, isoList = c(as.character(MBPIsoList$Isoform)), selfDf = MBPIsoList, gene = "NNAT")
otherDTE <- list()
for(t in WholeDESeqSig$age[WholeDESeqSig$age$isoform %in% redundant_MBP_iso,"isoform"]){
  otherDTE[[t]] <- plot_trans_exp_lifetime(t,class.files$glob_SQ,Exp$whole_group)
}
plot_grid(plotlist = otherDTE)

MBP <- list(
  tracks = ggTranPlots(inputgtf=input$merged_peptide_gtf, classfiles=proteinInput$t.class.files, isoList = c(as.character(MBPIsoList$Isoform)), selfDf = MBPIsoList, gene = "NNAT"),
  transcript = plot_trans_exp_lifetime("ONT18_5132_2313",class.files$glob_SQ,Exp$whole_group),
  protein = plot_trans_exp_lifetime("ONT18_5132_10030",class.files$glob_SQ,proteinCollapse$pWald$norm_counts)
)
plot_grid(MBP$tracks,
          plot_grid(
            plot_grid(MBP$transcript,MBP$protein,nrow=1,labels=c("B","C"),label_size = 25,scale = 0.9),
            plot_grid(plotlist = otherDTE[-1], nrow=2, labels = c("D"), label_size = 25,scale = 0.9),
            ncol=1, rel_heights = c(0.4,0.6)),scale=0.9)
plot_grid(MBP$tracks,plot_grid(MBP$transcript,MBP$protein,nrow=1,labels=c("B","C"),label_size = 25,scale = 0.9))

#merge(proteinCollapse$pWald$norm_counts[proteinCollapse$pWald$norm_counts$isoform == "ONT20_1826_53630",c("sample","normalised_counts")],
#      Exp$whole_group[Exp$whole_group$isoform == "ONT20_1826_53630",c("sample","normalised_counts")], by = "sample")

otherDTE <- list()
for(t in WholeDESeqSig$age[WholeDESeqSig$age$isoform %in% redundant_iso,"isoform"]){
  otherDTE[[t]] <- plot_trans_exp_lifetime(t,class.files$glob_SQ,Exp$whole_group)
}
plot_grid(plotlist = otherDTE)

pNNAT <- list(
  tracks = ggTranPlots(inputgtf=input$merged_peptide_gtf, classfiles=proteinInput$t.class.files, isoList = c(as.character(NNATIsoList$Isoform)), selfDf = NNATIsoList, gene = "NNAT"),
  transcript = plot_trans_exp_lifetime("ONT20_1826_53630",class.files$glob_SQ,Exp$whole_group),
  protein = plot_trans_exp_lifetime("ONT20_1826_53630",class.files$glob_SQ,proteinCollapse$pWald$norm_counts)
)


pdf(paste0(output_dir,"/NNAT.pdf"), width = 15, height = 10)
plot_grid(plot_grid(pNNAT$tracks,pNNAT$transcript,pNNAT$protein,nrow=1,rel_widths = c(0.5,0.25,0.25)), plot_grid(plotlist = otherDTE[-1],nrow=1),ncol=1, rel_heights = c(0.7,0.4))
dev.off()

plot_grid(pNNAT$tracks,
          plot_grid(
            plot_grid(pNNAT$transcript,pNNAT$protein,nrow=1,labels=c("B","C"),label_size = 25,scale = 0.9),
            plot_grid(plotlist = otherDTE[-1], nrow=2, labels = c("D"), label_size = 25,scale = 0.9),
            ncol=1, rel_heights = c(0.4,0.6)),scale=0.9)


