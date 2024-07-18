#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for main figures for isoform developmental paper
##
## ---------------------------------
## fig1:
## fig2: 
## fig3:
## fig4: Targeted dataset


SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir = paste0(SC_ROOT,"/Paper_Figures/outputFigs")


## ------ Figure 1: Whole dataset descriptives ------ 

fig1Whole <- list(
  # number of isoforms 
  #numIso = total_num_iso(class.files$glob_SQ, input_title = "", glimit = 10),
  # number of genes
  # class.files$glob_SQ %>% group_by(associated_gene) %>% tally()
  #numIsogene = numIsoGene(class.files$glob_SQ %>% filter(novelGene == "Annotated Genes")),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$glob_SQ_annoGene,rotate=TRUE) + theme(axis.title.y=element_blank(),
                                                                             axis.text.y=element_blank(),
                                                                             axis.ticks.y=element_blank())
  # cpat
  #cpat = plot_cpat(Cpat$whole, class.files$glob_SQ)
)



## ------ Figure 2: Targeted dataset descriptives ------ 

fig2Targeted <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$targ_SQ, input_title = "", glimit = 10),
  # number of genes
  numIsogene = numIsoGene(class.files$targ_SQ),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$targ_SQ) 
)


## ------ Figure 3: Differential expression analysis ------ 
tCol <- wes_palette("Royal1")[4]
fig3Diff <- list(
  # volcano plot across development
  volGroup = plot_volcano(diff_results=WholeDTE$age)[[1]],
  # top-ranked box-plot
  topRankedGroup = plot_trans_exp_individual("ONT18_5132_2319",class.files$glob_SQ,Exp$whole_group,"group", sqrt=TRUE, colourdots = alpha("#F8766D",0.3)),
  topRankedGroup2 = plot_trans_exp_individual("ONT18_5132_2313",class.files$glob_SQ,Exp$whole_group,"group", sqrt=TRUE, colourdots = "#00BFC4"),
  ##structuralcategory: class.files$glob_SQ[class.files$glob_SQ$isoform == "ONT18_5132_2319",]
  topRankedGroupVis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
              isoList = c("ONT18_5132_2319",RefIsoforms$MBP$transcript_id),
              colours = c(alpha("#F8766D",0.3),"black"), 
              simple=TRUE),
  topRankedGroupVis2 = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                  isoList = c("ONT18_5132_2313",RefIsoforms$MBP$transcript_id),
                                  colours = c("#00BFC4","black"), 
                                  simple=TRUE),
  volSex = plot_volcano(diff_results=WholeDTEAll$sex,interaction = "sex", chromosome = "allosomal")[[1]]
)

Add3Exp <- plot_trans_exp_individual("ONT10_4920_1919",class.files$glob_SQ,Exp$whole_group,"sex", sqrt=TRUE, colourdots = "#00BFC4")
##structuralcategory: class.files$glob_SQ[class.files$glob_SQ$isoform == "ONT18_5132_2319",]
Add3Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                isoList = c("ONT10_4920_1919",RefIsoforms$ADD3$transcript_id),
                                colours = c("#00BFC4","black"), 
                                simple=TRUE)

fig3Diff <- list(
  # volcano plot across development
  volGroup = plot_volcano(diff_results=WholeDTE$age)[[1]],
  # volcano plot across sex
  volSex = plot_volcano(diff_results=WholeDESeq$sex, interaction = "sex")[[1]],
  volSexSupp = plot_volcano(diff_results=WholeDESeq$sex, interaction = "sex", chromosome = "allosomal")[[1]],
  # top-ranked scatter plot
  topRankedGroup = plot_trans_exp_lifetime("ONT18_5132_2313",class.files$glob_SQ,Exp$whole_group),
  topRankedGroupVis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                   isoList = c("ONT18_5132_2313",RefIsoforms$MBP$transcript_id),
                                   colours = c(wes_palette("Royal1")[4],"black"), 
                                   simple=TRUE),
  topRankedGroup3 = plot_trans_exp_lifetime("ONT18_5132_2319",class.files$glob_SQ,Exp$whole_group),
  topRankedGroupVis3 = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                  isoList = c("ONT18_5132_2319",RefIsoforms$MBP$transcript_id),
                                  colours = c(wes_palette("Royal1")[4],"black"), 
                                  simple=TRUE),
  # top ranked antisense scatter plot 
  topRankedGroup2 = plot_trans_exp_lifetime("ONT8_3512_1715",class.files$glob_SQ,Exp$whole_group),
  topRankedGroup2Vis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                  isoList = c("ONT8_3512_1715",RefIsoforms$VXN$transcript_id),
                                  colours = c(wes_palette("Royal1")[4],"black"), 
                                  simple=TRUE),
  
  # top-ranked box-plot sex autosomal and allosomal chromosome
  #View(WholeDESeqSig$sex)
  #View(WholeDESeqSig$sex %>% filter(!chrom %in% c("chrY","chrX")))
  topRankedSex1 = plot_trans_exp_individual("ONTY_67_2",class.files$glob_SQ,Exp$whole_sex,"sex") + scale_fill_manual(values = c(tCol,tCol)),
  topRankedSex2 = plot_trans_exp_individual("ONTX_3476_23182",class.files$glob_SQ,Exp$whole_sex,"sex") + scale_fill_manual(values = c(tCol,tCol)),
  topRankedSex3 = plot_trans_exp_individual("ONT10_4920_1919",class.files$glob_SQ,Exp$whole_sex,"sex") + scale_fill_manual(values = c(tCol,tCol)),
  topRankedSex4 = plot_trans_exp_individual("ONT5_1265_990",class.files$glob_SQ,Exp$whole_sex,"sex") + scale_fill_manual(values = c(tCol,tCol)), 
  
  topRankedSex1Vis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
              isoList = c("ONTY_67_2",RefIsoforms$RPS4Y1$transcript_id),
              colours = c(wes_palette("Royal1")[4],"black"),
              simple=TRUE), 
  topRankedSex2Vis  = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
              isoList = c("ONTX_3476_23182",RefIsoforms$XIST$transcript_id),
              colours = c(wes_palette("Royal1")[4],"black"),
              simple=TRUE), 
  topRankedSex3Vis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
              isoList = c("ONT10_4920_1919",RefIsoforms$ADD3$transcript_id),
              colours = c(wes_palette("Royal1")[4],"black"),
              simple=TRUE),
  topRankedSex4Vis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                                 isoList = c("ONT5_1265_990",RefIsoforms$XIST$transcript_id),
                                 colours = c(wes_palette("Royal1")[4],"black"),
                                 simple=TRUE) 
)


## ------ Figure 4: GRIA3 ------ 

Fig4Targeted <- list(
  #volGroup = plot_volcano(diff_results=TargetedDESeqSig$age %>% filter(associated_gene %in% TargetGeneSZASD$targetGene))[[1]],
  topRankedUp = plot_trans_exp_lifetime("ONT7_6418_29474",class.files$glob_SQ,Exp$targeted_group),
  topRankedUpVis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                                 isoList = c("ONT7_6418_29474",RefIsoforms$CNTNAP2$transcript_id),
                                 colours = c(wes_palette("Cavalcanti1")[5],"#0C0C78"),
                                 simple=TRUE),
  topRankedDown = plot_trans_exp_lifetime("ONT18_407_1413",class.files$glob_SQ,Exp$targeted_group),
  topRankedDownVis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                               isoList = c("ONT18_407_1413",RefIsoforms$ANKRD12$transcript_id),
                               colours = c(wes_palette("Cavalcanti1")[5],"#0C0C78"),
                               simple=TRUE),
  GRIA3DIU = plotIFTargetedbyGene("GRIA3", paste0(dirnames$DIU,"targeted/group")), #+ 
    #scale_colour_manual(values = c(wes_palette("Royal1")[4],wes_palette("Royal2")[5])),
  GRIA3DGE = plot_trans_exp_lifetime(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$targeted_group,gene="GRIA3"),
  GRIA3DTE1 = plot_trans_exp_lifetime("ONTX_7115_9753",class.files$glob_targ_SQ,Exp$targeted_group),
  GRIA3DTE2 = plot_trans_exp_lifetime("ONTX_7115_973",class.files$glob_targ_SQ,Exp$targeted_group)
)



## ------ DTE across sex and development ------ 

WholeDTEsexAge <- lapply(intersect(WholeDTE$sex$isoform, WholeDTE$age$isoform), function(x) plot_trans_exp_individual(x,class.files$glob_SQ,Exp$whole_group,"both", sqrt=TRUE))


## ------ Output ------ 

# Figure 1
plot_grid(fig1Whole$numIso,fig2Targeted$numIso, labels = c("A","B"))
plot_grid(fig1Whole$numIsogene,fig2Targeted$numIsogene, labels = c("A","B"))
plot_grid(fig1Whole$numIsoCate,fig2Targeted$numIsoCate, labels = c("A","B"))

Fig3 <- plot_grid( 
          plot_grid(fig3Diff$volGroup, 
                    plot_grid(fig3Diff$topRankedGroupVis, fig3Diff$topRankedGroup,rel_heights = c(0.2,0.8),ncol=1), nrow = 1, rel_widths = c(0.4,0.6),labels = c("A","B")),
          plot_grid(fig3Diff$volSex,
                    plot_grid(fig3Diff$topRankedSex1Vis, fig3Diff$topRankedSex1,rel_heights = c(0.2,0.8),ncol=1),
                    plot_grid(fig3Diff$topRankedSex2Vis, fig3Diff$topRankedSex2,rel_heights = c(0.2,0.8),ncol=1),
                    plot_grid(fig3Diff$topRankedSex3Vis, fig3Diff$topRankedSex3,rel_heights = c(0.2,0.8),ncol=1),
                    plot_grid(fig3Diff$topRankedSex3Vis, fig3Diff$topRankedSex4,rel_heights = c(0.2,0.8),ncol=1), nrow = 1, rel_widths = c(0.4,0.2,0.2,0.2), labels = c("C","D","E","F")),
          ncol = 1
)

plot_grid(
  plot_grid(Fig4Targeted$volGroup,
            plot_grid(Fig4Targeted$topRankedUpVis, Fig4Targeted$topRankedUp, rel_heights = c(0.2,0.8),ncol=1),
            plot_grid(Fig4Targeted$topRankedDownVis, Fig4Targeted$topRankedDown, rel_heights = c(0.2,0.8),ncol=1), nrow = 1),
  plot_grid(Fig4Targeted$GRIA3DIU, Fig4Targeted$GRIA3DGE, Fig4Targeted$GRIA3DTE1, Fig4Targeted$GRIA3DTE2),
  ncol = 1
)

# targeted dataset
Top10TargetedDESeq2Sex = plot_top_results(TargetedDESeq2$sex,Exp$targeted,"sex")
Top10TargetedDESeq2Age = plot_top_results(TargetedDESeq2$age,Exp$targeted,"group")

TopRankedTargetedDESeq2Age <- list(
  box1 = plot_trans_exp_individual("ONT17_1642_1824",class.files$glob_targ_SQ,Exp$targeted_group,"group"), 
  box2 = plot_trans_exp_individual("ONT17_1642_18271",class.files$glob_targ_SQ,Exp$targeted_group,"group"),
  scatter1 = plot_trans_exp_lifetime("ONT17_1642_1824",class.files$glob_targ_SQ,Exp$targeted_group),
  scatter2 = plot_trans_exp_lifetime("ONT17_1642_18271",class.files$glob_targ_SQ,Exp$targeted_group)
)
TopRankedTargetedDESeq2AgeTracks <- list(
  RTN4 = ggTranPlots(gtf$merged,class.files$glob_targ_SQ,
                        isoList = c("ONT2_3232_1551",c("ENST00000394609.6","ENST00000405240.5","ENST00000357732.8")),
                        colours = c(wes_palette("Cavalcanti1")[5],rep("#0C0C78",length(RefIsoforms$RTN4)+1)), 
                        lines =  c(wes_palette("Cavalcanti1")[5],rep("#0C0C78",length(RefIsoforms$RTN4)+1)), 
                        gene = "SEPTIN4",simple=TRUE),
  
  SEPTIN4 = ggTranPlots(gtf$merged,class.files$glob_targ_SQ,
                                                  isoList = c("ONT17_1642_18247","ONT17_1642_18271",c("ENST00000583114.5","ENST00000672673.2","ENST00000580791.1")),
                                                  colours = c(wes_palette("Darjeeling2")[2],wes_palette("Darjeeling2")[3],rep("#0C0C78",length(RefIsoforms$SEPTIN4)+1)), 
                                                  lines =  c(wes_palette("Cavalcanti1")[5],wes_palette("Cavalcanti1")[4],rep("#0C0C78",length(RefIsoforms$SEPTIN4)+1)), 
                                                  gene = "SEPTIN4",simple=TRUE)
  
  
)




# whole transcriptome
Top10WholeDESeq2Sex = plot_top_results(WholeDESeq2$sex,Exp$wholeanno,"sex")
Top10WholeDESeq2Age = plot_top_results(WholeDESeq2$age,Exp$wholeanno,"group")


message("Number of significant differentially expressed transcripts across post and pre-natal: ", nrow(WholeDESeq2$age %>% filter(padj < 0.05)))
message("Number of significant differentially expressed transcripts across sex: ", nrow(WholeDESeq2$sex %>% filter(padj < 0.05)))



## ------ Figure 4: CPAT, alternative splicinge events ------ 

# tracks
ggTranPlots(inputgtf=gtf$glob_targ,classfiles=class.files$glob_targ_SQ,isoList=c("ONT1_2_13","ONT1_2_41"),colours = c("red","blue"),simple=TRUE)


## ------ Output ------ 

pdf(paste0(output_dir,"/WholeDescriptive.pdf"), width = 10, height = 10)
fig1Whole$numIsogene
fig1Whole$numIso
fig1Whole$numIsoCate
dev.off()

pdf(paste0(output_dir,"/TargetedDescriptive.pdf"), width = 10, height = 10)
fig2Targeted$numIsogene
fig2Targeted$numIso
fig2Targeted$numIsoCate
dev.off()

pdf("DTESexDevelopment_Supp.pdf", width = 12, height=8)
plot_grid(plotlist = WholeDTEsexAge)
dev.off()

pdf("DevelopmentVolcano.pdf", width = 8, height=6)
fig3Diff$volGroup
dev.off()

pdf("SexVolcano.pdf", width = 8, height=6)
fig3Diff$volSex
dev.off()

pdf("DevelopmentTopRanked.pdf", width = 6, height=6)
plot_grid(fig3Diff$topRankedGroupVis,fig3Diff$topRankedGroup,ncol=1, rel_heights = c(0.4,0.6))
dev.off()

pdf("DevelopmentTopRanked2.pdf", width = 6, height=6)
plot_grid(fig3Diff$topRankedGroupVis2, fig3Diff$topRankedGroup2,ncol=1, rel_heights = c(0.4,0.6))
dev.off()

pdf("SexTopRanked.pdf", width = 6, height=6)
plot_grid(Add3Vis,Add3Exp,ncol=1, rel_heights = c(0.4,0.6))
dev.off()

pdf(paste0(output_dir,"/targeted_DTE_group.pdf"), width = 20, height = 10)
plot_grid(TopRankedTargetedDESeq2Age$scatter1,TopRankedTargetedDESeq2AgeTracks$RTN4)
plot_grid(TopRankedTargetedDESeq2Age$scatter2,TopRankedTargetedDESeq2AgeTracks$SEPTIN4)
dev.off()


pdf(paste0(output_dir,"/Differential.pdf"), width = 20, height = 20)
Fig3
dev.off()

plot_volcano(diff_results=WholeDESeq$sex,interaction="sex")

p <- plot_volcano(diff_results=WholeDESeqSig$age)
ggsave(file="WholeDeSeq2Age.png", dpi=400)

plot_volcano(diff_results=TargetedDESeqSig$sex)
ggsave(file="TargetedDeSeq2Sex.png", dpi=400)

plotTargetedAge <- plot_volcano(diff_results=TargetedDESeq2$age)
plot_grid(plotTargetedAge[[1]])


DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% mutate(DGE_Dev = ifelse(Gene %in% WholeDESeqGeneSig$age$associated_gene,TRUE,FALSE),
                                                    DGE_Sex = ifelse(Gene %in% WholeDESeqGeneSig$sex$associated_gene,TRUE,FALSE))

head(DIUSig$wholeAllAge %>% filter(DGE_Dev == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, -as.numeric(totalChange)))
message("Number of genes with significant DIU across development: ", nrow(DIUSig$wholeAllAge))
message("Number of genes with significant DIU across development and podium Change: ", nrow(DIUSig$wholeAllAge %>% filter(podiumChange == TRUE)))
message("Number of genes with significant DIU across development and not podium Change: ", nrow(DIUSig$wholeAllAge %>% filter(podiumChange == FALSE)))
message("Number of genes with significant DIU across development and not DGE: ", nrow(DIUSig$wholeAllAge %>% filter(DGE_Dev == FALSE)))


DIU <- na.omit(DIUSig$wholeAllAge$Gene)
DGE <- na.omit(DIUSig$wholeAllAge[DIUSig$wholeAllAge$DGE_Dev == TRUE,"Gene"])
Podium <- na.omit(DIUSig$wholeAllAge[DIUSig$wholeAllAge$podiumChange == TRUE,"Gene"])


suppressMessages(library(VennDiagram))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

threevenndiagrams <- function(set1, set2,set3, name1, name2, name3){
  p <- venn.diagram(x = list(set1,set2, set3), 
                    label_alpha = 0, category.names = c(name1,name2, name3),filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', 
                    fill = c("#B3E2CD", "#FDCDAC","red"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                    print.mode = "raw")
  return(p)
}

p <- threevenndiagrams(DIU,DGE,Podium,"DIU","DGE","Major Isoform Switch")
plot_grid(grobTree(p))



plot_dendro_Tgene(dirnames$ficle,"DLGAP5")
DLGAP5ES <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"ES_events_counts.csv") %>% arrange(-numEvents) 
DLGAP5NE <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"NE_transcript_counts.csv") %>% arrange(-numNovelExons) 
DLGAP5IR <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"),"IR_events_counts.csv") %>% arrange(-numEvents) 
DLGAP5A5A3 <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/DLGAP5"), "DLGAP5_A5A3_transcript_counts.csv") %>% arrange(-numEvents)
DLGAP5All <- read.csv(paste0(dirnames$ficle,"/DLGAP5/Stats/", "DLGAP5_final_transcript_classifications.csv"))
DLGAP5AF <- DLGAP5All[DLGAP5All$AF == 1,"isoform"]
DLGAP5AT <- DLGAP5All[DLGAP5All$AT == 1,"isoform"]


table(FICLE_class_final_exp)

# DLGAP5
DLGAP5Prenatal <- class.files$glob_SQ[class.files$glob_SQ$associated_gene == "DLGAP5" & class.files$glob_SQ$DevStatus == "prenatal",] %>% arrange(-preReads)
DLGAP5Prenatal[DLGAP5Prenatal$isoform %in%  c("ONT14_1685_3519","ONT14_1685_3465","ONT14_1685_3463", "ONT14_1685_3452", "ONT14_1685_3439"),]
intersect(DLGAP5Prenatal$isoform,DLGAP5A5A3$transcriptID)
DLGAP5Iso <- data.frame(
  Isoform = unlist(DLGAP5Iso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "DLGAP5","transcript_id"])),
    `Pre-ES` = as.character(DLGAP5ES[,c("transcriptID")][1:4]),
    `Pre-NE` = as.character(DLGAP5NE[,c("transcriptID")][2:5]),
    `Pre-IR` = as.character(DLGAP5IR$transcriptID),
    `Pre-AF` = c("ONT14_1685_3519","ONT14_1685_3255","ONT14_1685_3463", "ONT14_1685_3439", "ONT14_1685_3180"),
    `Pre-AT` = c("ONT14_1685_3369", "ONT14_1685_3371", "ONT14_1685_3465","ONT14_1685_3273", "ONT14_1685_3190")
  )),
  Category = rep(names(DLGAP5Iso), lengths(DLGAP5Iso))
)
DLGAP5Iso$colour <- NA

DLGAP5ggTransFig <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
            isoList = c(as.character(DLGAP5Iso$Isoform)),
            selfDf = DLGAP5Iso, gene = "DLGAP5")

pdf("DLGAP5.pdf", width = 10, height = 6)
DLGAP5ggTransFig
dev.off()


dat <- DIUSig$wholeAllAge %>% mutate(totalChange = as.numeric(totalChange)) %>% filter(DGE_Dev == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, as.numeric(totalChange))
MORF4L2 <- list(
  IF = plotIFWholebyGene("MORF4L2",dirnames$DIU,facetTranscripts=TRUE)[[1]],
  #GeneExp = plot_trans_exp_individual(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGenes$whole_group,gene="MORF4L2",var="group") +
  #  scale_fill_manual(values = c(label_colour("Prenatal"),label_colour("Postnatal"))) + labs(title = NULL) ,
  #T1Exp = plot_trans_exp_lifetime("ONTX_6127_19264",class.files$glob_targ_SQ,Exp$whole_group,colpoints = wes_palette("Darjeeling1")[2]),
  #T2Exp = plot_trans_exp_lifetime("ONTX_6127_19265",class.files$glob_targ_SQ,Exp$whole_group,colpoints = wes_palette("Darjeeling2")[3]),
  Track = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                      isoList = c("ONTX_6127_19273","ONTX_6127_19269","ONTX_6127_19267","ONTX_6127_19265","ONTX_6127_19264",RefIsoforms$MORF4L2$transcript_id),
                      colours = c(rep("#00BFC4",2),alpha("#F8766D",0.5),rep("#00BFC4",2),"black"),
                      simple=TRUE)
)

DIUSummarySex <- DIUSig$wholeAllSex %>% mutate(totalChange = as.numeric(totalChange)) %>% filter(DGE_Sex == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, as.numeric(totalChange))
DIUSummarySex[DIUSummarySex$DGE_Dev == FALSE & DIUSummarySex$DGE_Sex == FALSE,]

GNAS_IF <- list(
  IF = plotIFWholebyGene("GNAS",dirnames$DIU,facetTranscriptsFeature=TRUE,sexFeature=TRUE)[[1]],
  vis = ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
              isoList = rev(c("ENST00000371085.8","ONT20_3069_10858","ONT20_3069_7102","ONT20_3069_9130","ONT20_3069_11037","ONT20_3069_10787")),
              colours = rev(c("black",alpha("#00BFC4",0.3),rep("#00BFC4",2),rep(alpha("#00BFC4",0.3),2))), 
              simple=TRUE)
)

pdf("SexDIU.pdf", width = 12, height = 10)
plot_grid(GNAS_IF$vis, GNAS_IF$IF, ncol = 1, rel_heights = c(0.4,0.6))
dev.off()

pdf("DevelopmentDIU.pdf", width = 12, height = 10)
plot_grid(MORF4L2$Track, MORF4L2$IF, ncol = 1, rel_heights = c(0.4,0.6))
dev.off()

