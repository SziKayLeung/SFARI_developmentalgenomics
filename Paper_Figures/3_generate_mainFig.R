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

suppressMessages(library("ggplot2"))
suppressMessages(library("wesanderson"))

SC_ROOT= "C:/Users/sl693/Dropbox/Scripts/SFARI_developmentalgenomics"
SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics"
source(paste0(SC_ROOT,"/Paper_Figures/SFARI_config.R"))
source(paste0(SC_ROOT,"/Paper_Figures/0_source_functions.R"))
output_dir <- "C:/Users/sl693/Dropbox/Scripts/SFARI_developmentalgenomics/output/"

figPlots <- list()

## ------ Figure 1: Whole dataset descriptives ------ 

pNum <- summary_plots()

# prenatal vs postnatal 
#geneNum <- read.csv(paste0(output_dir,"NumberofTranscriptsPrevsPost.csv"))
#ggplot(geneNum, aes(x = prenatalTranscripts, y = postnatalTranscripts)) + geom_point() +
#  theme_classic() + labs(x = "Number of transcripts: prenatal", y = "Number of transcripts: postnatal")

figPlots$structuralCategorySplit <- class.files$glob_SQ_annoGene %>% group_by(structural_category, DevStatus) %>% tally() %>% 
  mutate(DevStatus = factor(DevStatus,levels = c("prenatal","postnatal","Both"))) %>%
  ggplot(., aes(x = structural_category, y = n, fill = DevStatus)) + 
  geom_bar(stat="identity", position = position_dodge()) +
  mytheme + scale_fill_manual(name = "", 
                              labels = c("Prenatal","Postnatal","Both"), 
                              values = c(wes_palette("Royal1")[2],wes_palette("Royal1")[1],wes_palette("Royal1")[4])) +
  labs(x = "Structural category", y = "Number of transcripts of annotated genic features") +
  theme(legend.position = "top")


# comparison of whole vs targeted datasets across matched samples
figPlots$comp = whole_vs_targeted_plots(classfiles=class.files$glob_targ_SQ, 
                               wholeSamples=wholematchedsamples, 
                               targetedSamples=manifest[manifest$ID %in% targetedmatchedsamples,"Sample"], 
                               targetGene=selectedTargetGenes)
figPlots$comp + theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

# note there will be NA reads due to isoforms that are detected in the other samples that are not matching in the whole+targeted
# unique(matchedSumTargeted[matchedSumTargeted$dataset == "NA",])

# novel vs known transcripts: difference in expression
novelMean <- annoGenesStats$novelTrans %>% select(contains("Whole", ignore.case = FALSE)) %>% apply(., 1, mean) 
knownMean <- annoGenesStats$annoTrans %>% select(contains("Whole", ignore.case = FALSE)) %>% apply(., 1, mean) 
dat <- rbind(reshape2::melt(novelMean) %>% mutate(associated_transcript = "novel"), reshape2::melt(knownMean) %>% mutate(associated_transcript = "known"))
figPlots$novelKnownTranscriptExp <- ggplot(dat, aes(x = value, fill = associated_transcript)) + 
  geom_histogram(alpha = 0.3, position = "identity") + 
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "FL mean read count across all samples", y  = "Frequency", fill = "Transcript") +
  scale_x_continuous(trans='log10')


## ------ Figure 2: Targeted dataset descriptives ------ 

fig2Targeted <- list(
  # number of isoforms 
  numIso = total_num_iso(class.files$targ_SQ, input_title = "", glimit = 10),
  # number of genes
  numIsogene = numIsoGene(class.files$targ_SQ),
  # number of isoforms by category
  numIsoCate = plot_structural_cate(class.files$targ_SQ) 
)

figPlots$DLGAP5Dendro <- plot_dendro_Tgene(paste0(root_dir,"/ficle"),"DLGAP5")
DLGAP5Iso <- data.frame(
  Isoform = unlist(DLGAP5Iso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "DLGAP5", "transcript_id"])),
    `ES` = c("ONT14.1780.2914", "ONT14.1780.2873", "ONT14.1780.2793", "ONT14.1780.2671"),  # Updated ES IDs
    `NE` = c("ONT14.1780.3279", "ONT14.1780.3278", "ONT14.1780.2735", "ONT14.1780.2630"),  # Updated NE IDs
    `IR` = c("ONT14.1780.2615"),  # Updated IR ID
    `AF` = c("ONT14.1780.3457", "ONT14.1780.3403", "ONT14.1780.3380", "ONT14.1780.3193", "ONT14.1780.3117"),  # Updated AF IDs
    `AL` = c("ONT14.1780.3405", "ONT14.1780.3308", "ONT14.1780.3306", "ONT14.1780.3211", "ONT14.1780.3127")  # Updated AL IDs
  )),
  Category = rep(names(DLGAP5Iso), lengths(DLGAP5Iso))
)
DLGAP5Iso$colour <- NA
figPlots$DLGAP5Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ, isoList = c(as.character(DLGAP5Iso$Isoform)), 
                                  selfDf = DLGAP5Iso, gene = "DLGAP5")

RSP27A <- tabulateIF(class.files$glob_SQ %>% filter(associated_gene %in% "RPS27A"), "Whole") %>% 
  mutate(label = ifelse(perc > 10, isoform, NA), perc = perc / 100) %>% 
  plotIFGenes(.) + mytheme +
  labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme(legend.position = "right") + labs(y = "Isoform fracction")

RSP27AVis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                         isoList = c("ONT2.3331.23546","ONT2.3331.23589","ENST00000272317.11"),
                         colours = c(rep(alpha("#F8766D",0.3),2),"black"),
                         simple=TRUE)

## ------ Differential expression analysis ------

figPlots$Diff <- list()

# volcano plot across development
figPlots$Diff$volGroup = plot_volcano(diff_results=WholeDTE$age)[[1]]

# volcano plot across sex
figPlots$Diff$volSex = plot_volcano(diff_results=WholeDTEAll$sex, interaction = "sex")[[1]]
figPlots$Diff$volSexAllo = plot_volcano(diff_results=WholeDTEAll$sex, interaction = "sex", chromosome = "allosomal")[[1]]


# differential expressed across age
figPlots$Diff$topRankedGroup <- plot_trans_exp_individual("ONT18.5258.1932",class.files$glob_SQ,Exp$whole_group,"group", 
                                                          sqrt=TRUE, colourdots = alpha("#00BFC4",0.3))

figPlots$Diff$topRankedGroupVis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ, 
                                               isoList = c("ONT18.5258.1932",RefIsoforms$MBP$transcript_id), 
                                               colours = c("#00BFC4","black"), simple = TRUE)

figPlots$Diff$topRankedGroup2 <- plot_trans_exp_individual("ONT12.2697.32229",class.files$glob_SQ,Exp$whole_group,"group", 
                                                          sqrt=TRUE, colourdots = alpha("#F8766D",0.3))

figPlots$Diff$topRankedGroup2Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ, 
                                               isoList = c("ONT12.2697.32229",RefIsoforms$CHN1$transcript_id), 
                                               colours = c(alpha("#F8766D",0.3),"black"), simple = TRUE)

# differential expressed by sex
figPlots$Diff$topRankedSex <- plot_trans_exp_individual("ONT10.5139.1910",class.files$glob_SQ,Exp$whole_group,"sex", 
                                                        sqrt=TRUE, colourdots = "#00BFC4")

figPlots$Diff$topRankedSexVis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
                                             isoList = c("ONT10.5139.1910",RefIsoforms$ADD3$transcript_id),
                                             colours = c("#00BFC4","black"), simple=TRUE)

plot_grid(figPlots$Diff$topRankedGroupVis, figPlots$Diff$topRankedGroup, ncol = 1, rel_heights = c(0.3,0.7))
plot_grid(figPlots$Diff$topRankedSexVis, figPlots$Diff$topRankedSex, ncol = 1, rel_heights = c(0.3,0.7))

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

figPlots$WholeDTEsexAge <- lapply(intersect(WholeDTE$sex$isoform, WholeDTE$age$isoform), function(x) plot_trans_exp_individual(x,class.files$glob_SQ,Exp$whole_group,"both", sqrt=TRUE))


## ------ Differential transcript usage ------ 

# Overall
DIU <- intersect(na.omit(DIUSig$wholeAge$Gene),protein_coding_genes$V1)
DGE <- intersect(na.omit(DIUSig$wholeAge[DIUSig$wholeAge$DGE_Dev == TRUE,"Gene"]), protein_coding_genes$V1)
Podium <- intersect(na.omit(DIUSig$wholeAge[DIUSig$wholeAge$podiumChange == TRUE,"Gene"]), protein_coding_genes$V1)
DTE <- intersect(unique(na.omit(WholeDTE$age$associated_gene)), protein_coding_genes$V1)
p <- fourvenndiagrams(DIU,DGE,Podium,DTE,"DTU","DGE","Major Isoform Switch","DTE")
figPlots$DTEDGEDIUVenn <- plot_grid(grobTree(p))

# across development
DIUSig$wholeAge %>% mutate(totalChange = as.numeric(totalChange)) %>% filter(DGE_Dev == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, as.numeric(totalChange)) %>% head
figPlots$Diff$DIUDev <- plotIFWholebyGene("GPM6A",DIUnormExp$GMP6A_group,facetTranscriptsFeature=TRUE,sexFeature=FALSE)[[1]] + mytheme + theme(legend.position = "None")
figPlots$Diff$DIUDev_Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                                        isoList = c("ONT4.13313.18545","ONT4.13313.3469","ONT4.13313.18565","ONT4.13313.3461","ONT4.13313.3457",
                                                    RefIsoforms$GPM6A$transcript_id),
                                        colours = c(alpha("#00BFC4",0.3), "#00BFC4" ,rep(alpha("#00BFC4",0.3),3),"black"),
                                        simple=TRUE)
ExpGene$GPM6A_whole_group <- ExpGene$whole_group %>% filter(associated_gene == "GPM6A") %>% merge(., phenotype, by="sample")
figPlots$Diff$DIUDev_GeneExp <- 
  plot_trans_exp_individual(classfiles=class.files$glob_targ_SQ,Norm_transcount=ExpGene$GPM6A_whole_group,gene="GPM6A",var="group",
                                                          colourdots = "black") +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(color="black", size=2, alpha=0.9) + labs(y = "Normalized counts (log10)")

plot_grid(plot_grid(figPlots$Diff$DIUDev_Vis, labels = c("A")),
          plot_grid(figPlots$Diff$DIUDev, figPlots$Diff$DIUDev_GeneExp, nrow = 1, labels = c("B","C")), 
          ncol = 1, rel_widths = c(0.8,0.2))


# by sex
DIUSummarySex <- DIUSig$wholeSex %>% mutate(totalChange = as.numeric(totalChange)) %>% filter(DGE_Sex == FALSE, podiumChange == TRUE) %>% dplyr::arrange(FDR, as.numeric(totalChange))
DIUSummarySex[DIUSummarySex$DGE_Dev == FALSE & DIUSummarySex$DGE_Sex == FALSE,]

figPlots$Diff$DIUSex <- plotIFWholebyGene("GNAS",DIUnormExp$GNAS_sex,facetTranscriptsFeature=TRUE,sexFeature=TRUE)[[1]]
figPlots$Diff$DIUSex_Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_SQ,
            isoList = rev(c("ENST00000371085.8","ONT20.3125.11081","ONT20.3125.7140","ONT20.3125.9226","ONT20.3125.11276","ONT20.3125.10225")),
            colours = rev(c("black",alpha("#00BFC4",0.3),rep("#00BFC4",2),alpha("#00BFC4",0.3),"#00BFC4")), 
            simple=TRUE)
ExpGene$GNAS_whole_group <- ExpGene$whole_group %>% filter(associated_gene == "GNAS") %>% merge(., phenotype, by="sample")
figPlots$Diff$DIUSex_GeneExp <- plot_trans_exp_individual(transcript=NULL,
                                                          gene="GNAS",var="sex",sqrt=FALSE, 
                                                          colourdots = "black",
                                                          classfiles=class.files$glob_SQ,Norm_transcounts=ExpGene$GNAS_whole_group) +
                                scale_fill_manual(values = c("gray","gray"))

plot_grid(plot_grid(figPlots$Diff$DIUSex_Vis, labels = c("A")),
          plot_grid(figPlots$Diff$DIUSex[[1]], figPlots$Diff$DIUSex_GeneExp, nrow = 1, labels = c("B","C")), 
          ncol = 1, rel_widths = c(0.8,0.2))


## ------ AS events from FICLE ------ 

figPlots$ASEvents <- finalTranscriptClassificationGene2  %>% 
  filter(AS %in% c("A5A3", "AF", "AT", "ES", "IR", "NE_1st", "NE_Int", "NE_Last")) %>%
  mutate(DevStatus = factor(DevStatus, levels = c("prenatal","postnatal","Both"))) %>%
  group_by(AS, DevStatus) %>% tally(Frequency) %>% 
  ggplot(.,aes(x = reorder(AS, -n), y = n, fill = DevStatus)) + geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "AS events", y = "Frequency") %>% 
  scale_fill_manual(labels = c("Prenatal","Postnatal","Both"), name = "", values = c(wes_palette("Royal1")[2],wes_palette("Royal1")[1],wes_palette("Royal1")[4])) +
  labs(x = "AS events", y = "Frequency") +
  theme(legend.position = "top") +
  theme_classic()


## ------ Targeted data ------ 

SFARI = c('ADNP','AGAP1','ANK2','ANKRD11','ARID1B','ASH1L','ASPM','AUTS2','BCL11A','CACNA1C','CACNA1G','CACNA1H','CADPS','CADPS2','CD38','CDH13','CDH2','CDK13','CELF6','CHD2','CHD8','CLTCL1','CNOT1','CNTN4','CNTN6','CNTNAP2','CSMD1','CTNNA2','CYFIP1','DAGLA','DCC','DISC1','DLG4','DLGAP2','DMD','DPYD','DRD2','DRD3','DYRK1A','ELP4','EP300','FBXO40','FMR1','FOXP1','FOXP2','GABBR2','GPC6','GRIA3','GRIK2','GRIK3','GRIN1','GRIN2A','GRIN2B','H1-4','HERC1','IL1RAPL1','IMMP2L','ITPR1','KAT6B','KATNAL2','KCTD13','KDM6A','KDM6B','KIRREL3','LRBA','MAPT','MCM4','MCPH1','MECP2','MED13L','MEIS2','MET','NEGR1','NF1','NFIA','NFIB','NKX2-2','NLGN1','NLGN2','NLGN3','NLGN4X','NR3C2','NRXN1','NTNG1','NXPH1','OXTR','PARD3B','PAX6','PBX1','PCDH10','PCDH9','PHB','PJA1','POGZ','PRKN','PTEN','PTPRT','QRICH1','RBFOX1','RELN','RPL10','RUNX1T1','SCN1A','SCN2A','SCN8A','SET','SETD1A','SEZ6L2','SHANK3','SLC4A10','SLC9A9','SON','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','TSHZ3','UBE3A','USP7','VPS13B','WWOX','ZBTB16','ZBTB20','ZMYM2','ZNF18','ZNF804A')

class.files$glob_targ_SQ_counts_matrix <- class.files$glob_targ_SQ %>% select(contains(c("Whole","Targeted"), ignore.case = FALSE)) 
colnames(class.files$glob_targ_SQ_counts_matrix) <- paste(colnames(class.files$glob_targ_SQ_counts_matrix),"Reads",sep="_") 
class.files$glob_targ_SQ_counts_matrix <- cbind(class.files$glob_targ_SQ[,c("isoform","associated_gene","structural_category")],class.files$glob_targ_SQ_counts_matrix)
CombinedSFARI <- tabulateIF(class.files$glob_targ_SQ_counts_matrix %>% filter(associated_gene %in% SFARI), "Reads")
figPlots$SFARIDIU <- plotIFGenes(CombinedSFARI)

##--- DAGLA ---
figPlots$DAGLADendro <- plot_dendro_Tgene(paste0(root_dir,"/ficle"),"DAGLA", cfiles = class.files$glob_targ_SQ)
DAGLAIso <- data.frame(
  Isoform = unlist(DAGLAIso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "DAGLA", "transcript_id"])),
    `ES` = as.character(c("ONT11.3444.7446", "ONT11.3444.7443", "ONT11.3444.7427", "ONT11.3444.7422")),  
    `AF, ES` = as.character(c("ONT11.3444.7620")),  
    `AF` = as.character(c("ONT11.3444.7645","ONT11.3444.6615")),  
    `AL` = as.character("ONT11.3444.7914"),  
    `AF, AL` = as.character(c("ONT11.3444.7617"))  
  )),
  Category = rep(names(DAGLAIso), lengths(DAGLAIso))
)
DAGLAIso$colour <- NA
figPlots$DAGLAVis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                                isoList = c(as.character(DAGLAIso$Isoform)),
                                selfDf = DAGLAIso, gene = "DAGLA")

##--- FOXP2 ---
figPlots$FOXP2PDendro <- plot_dendro_Tgene(paste0(root_dir,"/ficle"),"FOXP2", cfiles = class.files$glob_targ_SQ)
FOXP2Iso <- data.frame(
  Isoform = unlist(FOXP2Iso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "FOXP2", "transcript_id"])),
    `AP` = as.character(c("ONT7.4758.853", "ONT7.4758.46", "ONT7.4758.1615", "ONT7.4758.188", "ONT7.4758.194", "ONT7.4758.987")),
    `AP, NE` = as.character(c("ONT7.4758.4670", "ONT7.4758.4597", "ONT7.4758.3486", "ONT7.4758.3496", "ONT7.4758.3122", "ONT7.4758.34", "ONT7.4758.3243")),
    `AF, NE` = as.character(c("ONT7.4758.2928", "ONT7.4758.2695", "ONT7.4758.2801", "ONT7.4758.2809")),
    `AF` = as.character(c("ONT7.4758.2661", "ONT7.4758.2859"))
  )),
  Category = rep(names(FOXP2Iso), lengths(FOXP2Iso))
)
FOXP2Iso$colour <- NA
figPlots$FOXP2Vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                        isoList = c(as.character(FOXP2Iso$Isoform)),
                        selfDf = FOXP2Iso, gene = "FOXP2")

##--- CACNA1G ---
figPlots$CACNA1GPDendro <- plot_dendro_Tgene(paste0(root_dir,"/ficle"),"CACNA1G", cfiles = class.files$glob_targ_SQ)
CACNA1GIso <- data.frame(
  Isoform = unlist(CACNA1GIso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "CACNA1G", "transcript_id"])),
    `ES` = c("ONT17.1148.1890", "ONT17.1148.1880", "ONT17.1148.1809", "ONT17.1148.1773", "ONT17.1148.1772"),
    `IR` = c("ONT17.1148.1527", "ONT17.1148.1519", "ONT17.1148.1508", "ONT17.1148.1503", "ONT17.1148.1489"),
    `IR, AF` = c("ONT17.1148.1394", "ONT17.1148.1383", "ONT17.1148.1329", "ONT17.1148.1373")
  )),
  Category = rep(names(CACNA1GIso), lengths(CACNA1GIso))
)
CACNA1GIso$colour <- NA
figPlots$CACNA1Gvis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                          isoList = c(as.character(CACNA1GIso$Isoform)),
                          selfDf = CACNA1GIso, gene = "CACNA1G")

##--- ERCC ---

ercc_usage <- plot_usage_persample(ercc.class.file, ercc.demux)
figPlots$ERCC <- plot_grid(plotlist = ercc_usage$plotsFSM, ncol = 1, labels = c("A","B","C"))
ercc_usage_average <- plot_average_usage_across_all_samples(ercc_usage$tables) 
mean(
  mean(ercc_usage$tables$BC01[ercc_usage$tables$BC01$structural_category == "FSM","nn"]),
  mean(ercc_usage$tables$BC04[ercc_usage$tables$BC04$structural_category == "FSM","nn"]),
  mean(ercc_usage$tables$BC06[ercc_usage$tables$BC06$structural_category == "FSM","nn"])
)

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

pdf(paste0(output_dir,"/numKnownGenicFeatures_1c.pdf"), width = 10, height = 10)
pNum[[2]]
dev.off()

pdf(paste0(output_dir,"/numKnownGenicFeatures_1d_1e.pdf"), width = 15, height = 10)
pNum[[1]]
dev.off()


pdf(paste0(output_dir,"/TargetedDescriptive.pdf"), width = 10, height = 10)
fig2Targeted$numIsogene
fig2Targeted$numIso
fig2Targeted$numIsoCate
dev.off()

pdf(paste0(output_dir,"DTESexDevelopment_Supp.pdf"), width = 14, height=8)
plot_grid(plotlist = figPlots$WholeDTEsexAge)
dev.off()

pdf(paste0(output_dir, "DevVolcano.pdf"), width = 12, height = 8)
figPlots$Diff$volGroup
dev.off()

pdf(paste0(output_dir, "SexVolcano.pdf"), width = 12, height = 8)
figPlots$Diff$volSex
figPlots$Diff$volSexAllo
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

pdf(paste0(output_dir,"/ERCC.pdf"), width = 10, height = 15)
figPlots$ERCC
dev.off()

plot_volcano(diff_results=WholeDESeq$sex,interaction="sex")

p <- plot_volcano(diff_results=WholeDESeqSig$age)
ggsave(file="WholeDeSeq2Age.png", dpi=400)

plot_volcano(diff_results=TargetedDESeqSig$sex)
ggsave(file="TargetedDeSeq2Sex.png", dpi=400)

plotTargetedAge <- plot_volcano(diff_results=TargetedDESeq2$age)
plot_grid(plotTargetedAge[[1]])

pdf("SexDIU.pdf", width = 12, height = 10)
plot_grid(GNAS_IF$vis, GNAS_IF$IF, ncol = 1, rel_heights = c(0.4,0.6))
dev.off()

# Figure 2
pdf(paste0(output_dir, "Figure2A_DLGAP5.pdf"), width = 10, height = 6)
figPlots$DLGAP5Vis
dev.off()

pdf(paste0(output_dir, "Figure2B_DLGAP5.pdf"), width = 5, height = 6)
figPlots$DLGAP5Dendro
dev.off()

# Figure 4
pdf(paste0(output_dir, "NovelTopRankedDevDiff.pdf"), width = 10, height= 8)
plot_grid(figPlots$Diff$topRankedGroup2Vis,figPlots$Diff$topRankedGroup2, ncol = 1, rel_widths = c(0.2,0.8))
dev.off()

pdf(paste0(output_dir,"Figure4C_GPM6A_DIU.pdf"),width = 10, height = 6)
figPlots$Diff$DIUDev_Vis
dev.off()

pdf(paste0(output_dir,"Figure4D_GPM6A_DIU.pdf"),width = 15, height = 6)
figPlots$Diff$DIUDev
dev.off()

pdf(paste0(output_dir,"Figure4E_GPM6A_DIU.pdf"),width = 5, height = 6)
figPlots$Diff$DIUDev_GeneExp
dev.off()

# Figure 5
pdf(paste0(output_dir, "Figure5A_SFARIDIU.pdf"), width = 30, height = 12)
figPlots$SFARIDIU + theme(legend.position = "bottom", text=element_text(size=20))
dev.off()

pdf(paste0(output_dir,"DAGLA_Figure5A.pdf"), width = 10, height = 6)
figPlots$DAGLAVis
dev.off()

pdf(paste0(output_dir,"FOXP2_Figure5B.pdf"), width = 10, height = 6)
figPlots$FOXP2Vis
dev.off()

pdf(paste0(output_dir,"CACNA1G_Figure5C.pdf"), width = 10, height = 6)
figPlots$CACNA1Gvis
dev.off()


# supplementary figures
pdf(paste0(output_dir,"FICLEASEventsPerGene.pdf"), width = 10, height = 6)
figPlots$ASEvents
dev.off()

pdf("RSP27A_Supplementary.pdf", width = 10, height = 6)
plot_grid(RSP27AVis,RSP27A, rel_widths = c(0.4,0.6))
dev.off()
