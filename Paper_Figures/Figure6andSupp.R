monoAllelicDDP = c('ADNP','ANK2','ANKRD11','ARID1B','ASH1L','AUTS2','BCL11A','BCL11B','CACNA1C','CACNA1G','CACNA1H','CDH1','CDH2','CDK13','CHD2','CHD8','CLCN3','CNOT1','DLG4','DMD','DYRK1A','EIF2S3','EP300','FMR1','FOXP1','FOXP2','GABBR2','GATA3','GATA4','GATA6','GLMN','GRIA3','GRIK2','GRIN1','GRIN2A','GRIN2B','H1-4','HNF1B','HNF4A','IL1RAPL1','ITPR1','KAT6B','KDM6A','KDM6B','KIRREL3','KMT2D','LMNA','MAGI2','MECP2','MED13L','MEIS2','MNX1','NF1','NFIA','NFIB','NLGN3','NLGN4X','NRXN1','OPA1','PAX6','PBX1','POGZ','POLA1','PTEN','QRICH1','RANBP2','RBFOX1','RPL10','SCN1A','SCN2A','SCN8A','SET','SETD1A','SHANK3','SHH','SLC9A9','SMARCE1','SON','SOX9','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','UBE3A','USP7','WFS1','ZBTB20','ZEB2','ZMYM2')

SFARI = c('ADNP','AGAP1','ANK2','ANKRD11','ARID1B','ASH1L','ASPM','AUTS2','BCL11A','CACNA1C','CACNA1G','CACNA1H','CADPS','CADPS2','CD38','CDH13','CDH2','CDK13','CELF6','CHD2','CHD8','CLTCL1','CNOT1','CNTN4','CNTN6','CNTNAP2','CSMD1','CTNNA2','CYFIP1','DAGLA','DCC','DISC1','DLG4','DLGAP2','DMD','DPYD','DRD2','DRD3','DYRK1A','ELP4','EP300','FBXO40','FMR1','FOXP1','FOXP2','GABBR2','GPC6','GRIA3','GRIK2','GRIK3','GRIN1','GRIN2A','GRIN2B','H1-4','HERC1','IL1RAPL1','IMMP2L','ITPR1','KAT6B','KATNAL2','KCTD13','KDM6A','KDM6B','KIRREL3','LRBA','MAPT','MCM4','MCPH1','MECP2','MED13L','MEIS2','MET','NEGR1','NF1','NFIA','NFIB','NKX2-2','NLGN1','NLGN2','NLGN3','NLGN4X','NR3C2','NRXN1','NTNG1','NXPH1','OXTR','PARD3B','PAX6','PBX1','PCDH10','PCDH9','PHB','PJA1','POGZ','PRKN','PTEN','PTPRT','QRICH1','RBFOX1','RELN','RPL10','RUNX1T1','SCN1A','SCN2A','SCN8A','SET','SETD1A','SEZ6L2','SHANK3','SLC4A10','SLC9A9','SON','SRRM2','STAG1','STXBP1','SYNGAP1','SYP','TBL1XR1','TBR1','TCF4','TRIO','TSC1','TSC2','TSHZ3','UBE3A','USP7','VPS13B','WWOX','ZBTB16','ZBTB20','ZMYM2','ZNF18','ZNF804A')

schemaGenes <- read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/schema_genes.txt", col.names = c("Gene"))
GWAS = c("ACTR1B", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CSMD1", "CTD-2008L17.2", "CUL9", "DCC", "DLGAP2", "DPYD", "EMX1", "ENOX1", "EPN2", "EYS", "FURIN", "GABBR2", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "IL1RAPL1", "IMMP2L", "IRF3", "KIAA1549", "KLF6", "LINC00320", "LINC01088", "LRRC4B", "MAD1L1", "MAN2A1", "MAPT", "MSI2", "NAB2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCGF3", "PCNXL3", "PDE4B", "PJA1", "PLCH2", "PTPRD", "R3HDM2", "RP11-399D6.2", "RP11-507B12.2", "SGCD", "SLC39A8", "SLC4A10", "SNAP91", "SP4", "THAP8", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WSCD2", "ZNF804A", "ZNF823", "ZNF835")

biallelicDDP = c('ASPM','CACNA1G','CHL1','CISD2','CLCN3','CNTNAP2','CTNNA2','DCC','EIF2AK3','EOMES','GLIS3','GPC6','GRIK2','GRIN1','GRIN2A','GRM1','HADH','HERC1','HPSE2','ITPR1','KDM6B','KIAA1109','LMNA','LRBA','MCPH1','NRXN1','NTNG1','ONECUT1','PDIA6','PDX1','PEX16','PPP1R15B','PTF1A','QARS1','RELN','RFX6','SLC39A8','SLF2','TRMT10A','TSPEAR','UBE3B','VPS13B','WFS1','WWOX','ZBTB16','ZFP57')

length(unique(as.character(monoAllelicDDP))) + length(unique(as.character(biallelicDDP)))
monoAllelicDDPUnique <- setdiff(monoAllelicDDP, SFARI)
schemaGenesAlone <- setdiff(setdiff(schemaGenes$Gene,SFARI),monoAllelicDDP)

length(unique(c(monoAllelicDDP, biallelicDDP, schemaGenes$Gene, SFARI)))


tabulateIF <- function(classf, countcol){
  
  Counts <- classf %>% select(isoform,contains(countcol))
  rownames(Counts) <- Counts$isoform
  Counts <- Counts %>% select(-isoform)
  
  
  # Calculate the mean of normalised expression across all the samples per isoform
  meandf <- data.frame(meanvalues = apply(Counts,1,mean)) %>%
    rownames_to_column("isoform") %>% 
    # annotate isoforms with associated_gene and structural category
    left_join(., classf[,c("isoform","associated_gene","structural_category")], by = "isoform")  
  
  # Group meandf by associated_gene and calculate the sum of mean values for each group
  grouped <- aggregate(meandf$meanvalues, by=list(associated_gene=meandf$associated_gene), FUN=sum)
  
  # Calculate the proportion by merging back, and divide the meanvalues by the grouped values (x)
  merged <- meandf %>% 
    left_join(grouped, by = "associated_gene") %>%
    mutate(perc = meanvalues / x * 100) 
  return(merged)
}

WholeMonoAllelicDDP <- tabulateIF(class.files$glob_SQ %>% filter(associated_gene %in% monoAllelicDDP ), "Whole")
TargetedMonoAllelicDDP <- tabulateIF(class.files$targ_SQ %>% filter(associated_gene %in% monoAllelicDDP), "Targeted")

# combined dataset
class.files$glob_targ_SQ_counts_matrix <- cbind(select(class.files$glob_targ_SQ_counts,contains("Whole")),select(class.files$glob_targ_SQ_counts,contains("Targeted")))
colnames(class.files$glob_targ_SQ_counts_matrix) <- paste(colnames(class.files$glob_targ_SQ_counts_matrix),"Reads",sep="_") 
class.files$glob_targ_SQ_counts_matrix <- cbind(class.files$glob_targ_SQ_counts[,c("isoform","associated_gene","structural_category")],class.files$glob_targ_SQ_counts_matrix)

CombinedSFARI <- tabulateIF(class.files$glob_targ_SQ_counts_matrix %>% filter(associated_gene %in% SFARI), "Reads")
CombinedMonoAllelicDDP <- tabulateIF(class.files$glob_targ_SQ_counts_matrix %>% filter(associated_gene %in% monoAllelicDDPUnique), "Reads")  
CombinedSchema <- tabulateIF(class.files$glob_targ_SQ_counts_matrix %>% filter(associated_gene %in% schemaGenesAlone), "Reads")  
      
plotIFGenes <- function(dat){
  dat <- dat %>% mutate(structural_category = factor(structural_category, levels = c("FSM","ISM","NIC","NNC", "Genic_Genomic")))
  p <- ggplot(dat, aes(x = associated_gene, y = as.numeric(perc), fill = forcats::fct_rev(structural_category))) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    #scale_color_manual(values = rep(NA, length(unique(minorgrouped$gene)))) + 
    labs(x = "Gene", y = "Isoform fraction (%)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  if(length(unique(dat$structural_category)) == 4){
    p <- p + scale_fill_manual(name = "Isoform Classification", values = rev(c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                               alpha("#F8766D",0.8),alpha("#F8766D",0.3)))) 
  }else{
    p <- p + scale_fill_manual(name = "Isoform Classification", values = rev(c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                               alpha("#F8766D",0.8),alpha("#F8766D",0.3),"gray"))) 
  }
    #+
    #theme(legend.position = "None")
  
  return(p)
}

library("cowplot")
plot_grid(plotIFGenes(CombinedSFARI), plotIFGenes(CombinedMonoAllelicDDP), plotIFGenes(CombinedSchema), rel_widths = c(0.7,0.2,0.10), nrow = 1)

RSP27A <- tabulateIF(class.files$glob_SQ%>% filter(associated_gene %in% "RPS27A"), "Whole") %>% 
  mutate(label = ifelse(perc > 10, isoform, NA), perc = perc / 100) %>% 
  plotIFGenes(.) + mytheme +
  labs(x = "") + scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme(legend.position = "right") + labs(y = "Isoform fracction")

RSP27AVis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
            isoList = c("ONT2_3223_23931","ONT2_3223_23974","ENST00000272317.11"),
            colours = c(rep(alpha("#F8766D",0.3),2),"black"),
            simple=TRUE)

pdf("RSP27A_Supplementary.pdf", width = 10, height = 6)
plot_grid(RSP27AVis,RSP27A, rel_widths = c(0.4,0.6))
dev.off()


IFCombined <- plotIFGenes(CombinedSFARI)
dirnames$ficleRerun <- "/lustre/home/sl693/FICLE/TargetGenes"
DAGLAp <- plot_dendro_Tgene("/lustre/home/sl693/FICLE/TargetGenes", "DAGLA", cfiles = class.files$glob_targ_SQ)
FOXP2p <- plot_dendro_Tgene("/lustre/home/sl693/FICLE/TargetGenes", "FOXP2", cfiles = class.files$glob_targ_SQ)
CACNA1Gp <- plot_dendro_Tgene("/lustre/home/sl693/FICLE/TargetGenes", "CACNA1G", cfiles = class.files$glob_targ_SQ)

pdf("Figure6Top.pdf", width = 30, height = 12)
#plot_grid(IFCombined + theme(legend.position = "None"), plot_grid(DAGLAp, FOXP2p, CACNA1Gp, nrow = 1), ncol = 1, rel_heights = c(0.7,0.3))
IFCombined + theme(legend.position = "None", text=element_text(size=20))
dev.off()

pdf("Figure6Mid.pdf", width = 30, height = 6)
plot_grid(DAGLAp, FOXP2p, CACNA1Gp, nrow = 1)
dev.off()

FOXP2ES <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/FOXP2"),"ES_events_counts.csv") %>% arrange(-numEvents) 
FOXP2NE <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/FOXP2"),"NE_transcript_counts.csv") %>% arrange(-numNovelExons) 
FOXP2All <- read.csv(paste0(dirnames$ficle,"/FOXP2/Stats/", "FOXP2_final_transcript_classifications.csv"))
FOXP2AF <- FOXP2All[FOXP2All$AF == 1,"isoform"]
FOXP2AllAT <-FOXP2All[FOXP2All$AT == 1,"isoform"]

FOXP2Iso <- data.frame(
  Isoform = unlist(FOXP2Iso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "FOXP2","transcript_id"])),
    `AP` = as.character(c("ONT7_4543_851", "ONT7_4543_46", "ONT7_4543_1616", "ONT7_4543_184", "ONT7_4543_190", "ONT7_4543_984" )),
    `AP, NE` = as.character(c("ONT7_4543_4529","ONT7_4543_4470","ONT7_4543_3377", "ONT7_4543_3385", "ONT7_4543_3061", "ONT7_4543_34", "ONT7_4543_3173")),
    `AF, NE` = as.character(c("ONT7_4543_2880", "ONT7_4543_2688", "ONT7_4543_2789", "ONT7_4543_2785")),
    `AF` = as.character(c("ONT7_4543_2660","ONT7_4543_2823"))
  )),
  Category = rep(names(FOXP2Iso), lengths(FOXP2Iso))
)
FOXP2Iso$colour <- NA
FOXP2vis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
            isoList = c(as.character(FOXP2Iso$Isoform)),
            selfDf = FOXP2Iso, gene = "FOXP2")



CACNA1GES <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/CACNA1G"),"ES_events_counts.csv") %>% arrange(-numEvents) 
CACNA1GNE <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/CACNA1G"),"NE_transcript_counts.csv") %>% arrange(-numNovelExons) 
CACNA1GIR <- input_FICLE_splicing_results(paste0(dirnames$ficle,"/CACNA1G"),"IR_transcript_level.csv") 
CACNA1GAll <- read.csv(paste0(dirnames$ficle,"/CACNA1G/Stats/", "CACNA1G_final_transcript_classifications.csv"))
CACNA1GAF <- CACNA1GAll[CACNA1GAll$AF == 1,"isoform"]
CACNA1GAllAT <-CACNA1GAll[CACNA1GAll$AT == 1,"isoform"]

CACNA1GIso <- data.frame(
  Isoform = unlist(CACNA1GIso <- list(
    Reference = as.character(unique(gtf$ref[gtf$ref$gene_id == "CACNA1G","transcript_id"])),
    #`SS` = as.character(class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene == "CACNA1G",] %>% arrange(length) %>% .[,c("isoform")])[1:5],
    `ES` = as.character(CACNA1GES$transcriptID[1:5]),
    `IR` = c("ONT17_1134_1521","ONT17_1134_1513","ONT17_1134_1502","ONT17_1134_1497","ONT17_1134_1483"),
    `IR AF` = c("ONT17_1134_1394","ONT17_1134_1384","ONT17_1134_1331","ONT17_1134_1375")
    
  )),
  Category = rep(names(CACNA1GIso), lengths(CACNA1GIso))
)
CACNA1GIso$colour <- NA
CACNA1Gvis <- ggTranPlots(inputgtf=gtf$merged,classfiles=class.files$glob_targ_SQ,
                          isoList = c(as.character(CACNA1GIso$Isoform)),
                          selfDf = CACNA1GIso, gene = "CACNA1G")

pdf("CACNA1G.pdf", width = 10, height = 6)
CACNA1Gvis
dev.off()

pdf("FOXP2.pdf", width = 10, height = 6)
FOXP2vis
dev.off()

pdf("DAGLA.pdf", width = 10, height = 6)
DAGLAvis
dev.off()

####
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% mutate(DTE_Dev = ifelse(Gene %in% WholeDTE$age$associated_gene, TRUE, FALSE),
                                                    IsoformSwitch = ifelse(podiumChange == TRUE, TRUE, FALSE))
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% filter(!is.na(Gene))

p <- venn.diagram(
  x = list(intersect(unique(WholeDTE$age$associated_gene), ProteinCodingGenes), 
           intersect(unique(WholeDESeqGeneSig$age$associated_gene), ProteinCodingGenes), 
           intersect(DIUSig$wholeAllAge$Gene, ProteinCodingGenes),
           intersect(DIUSig$wholeAllAge[DIUSig$wholeAllAge$podiumChange == TRUE, "Gene"], ProteinCodingGenes)),
  fill = c("red", "blue","green", "light green"),
  category.names = c("DTE" , "DGE" , "DIU", "Isoform Switch"),
  print.mode = "raw", filename = NULL, output = TRUE
)
plot_grid(p)


## final transcript FICLE
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% mutate(DTE_Dev = ifelse(Gene %in% WholeDTE$age$associated_gene, TRUE, FALSE),
                                                    IsoformSwitch = ifelse(podiumChange == TRUE, TRUE, FALSE))
DIUSig$wholeAllAge <- DIUSig$wholeAllAge %>% filter(!is.na(Gene))

p <- venn.diagram(
  x = list(intersect(unique(WholeDTE$age$associated_gene), ProteinCodingGenes), 
           intersect(unique(WholeDESeqGeneSig$age$associated_gene), ProteinCodingGenes), 
           intersect(DIUSig$wholeAllAge$Gene, ProteinCodingGenes),
           intersect(DIUSig$wholeAllAge[DIUSig$wholeAllAge$podiumChange == TRUE, "Gene"], ProteinCodingGenes)),
  fill = c("red", "blue","green", "light green"),
  category.names = c("DTE" , "DGE" , "DIU", "Isoform Switch"),
  print.mode = "raw", filename = NULL, output = TRUE
)
plot_grid(p)


finalTranscriptClassification <- distinct(fread(paste0(dirnames$ficle,"/all_final_transcript_classifications.csv"), data.table = F))
finalTranscriptClassification <- finalTranscriptClassification %>% filter(isoform %in% class.files$glob_SQ$isoform)
finalTranscriptClassification <- merge(finalTranscriptClassification, class.files$glob_SQ[,c("isoform", "associated_gene", "associated_transcript","structural_category")], by = "isoform")
finalTranscriptClassificationGene <- aggregate(. ~ associated_gene, finalTranscriptClassificationTranscript %>% select(-isoform), sum)



finalTranscriptClassificationGene  %>% filter(AS %in% c("A5A3", "AF", "AT", "ES", "IR", "NE_1st", "NE_Int", "NE_Last")) %>%
  group_by(AS) %>% tally(Frequency) %>% 
  ggplot(.,aes(x = reorder(AS, -n), y = n)) + geom_bar(stat = "identity") +
  mytheme + labs(x = "AS events", y = "Frequency") + facet_grid(~DevStatus)

finalTranscriptClassificationTranscript <- merge(finalTranscriptClassification, class.files$glob_SQ[,c("isoform","DevStatus")], by = "isoform")
finalTranscriptClassificationGene <- aggregate(. ~ associated_gene + DevStatus, 
                                               finalTranscriptClassificationTranscript %>% select(-isoform), 
                                               sum)
finalTranscriptClassificationGene2 <- reshape2::melt(finalTranscriptClassificationGene, variable.name = "AS", value.name = "Frequency", id = c("associated_gene","DevStatus"))


finalTranscriptClassificationGene2  %>% filter(AS %in% c("A5A3", "AF", "AT", "ES", "IR", "NE_1st", "NE_Int", "NE_Last")) %>%
  mutate(DevStatus = factor(DevStatus, levels = c("prenatal","postnatal","Both"))) %>%
  group_by(AS, DevStatus) %>% tally(Frequency) %>% 
  ggplot(.,aes(x = reorder(AS, -n), y = n, fill = DevStatus)) + geom_bar(stat = "identity", position = position_dodge()) +
  mytheme + labs(x = "AS events", y = "Frequency") %>% 
  scale_fill_manual(labels = c("Prenatal","Postnatal","Both"), name = "", values = c(wes_palette("Royal1")[2],wes_palette("Royal1")[1],wes_palette("Royal1")[4])) +
  labs(x = "AS events", y = "Frequency") +
  theme(legend.position = "top")

