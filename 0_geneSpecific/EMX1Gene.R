suppressMessages(library("rtracklayer"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))


source("/lustre/projects/Research_Project-MRC148213/sl693/scripts/LOGen/merge_characterise_dataset/run_ggtranscript.R")
dir <- "/lustre/projects/Research_Project-MRC148213/sl693/SFARI/"
gtf <- paste0(dir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf")
gtf <- rtracklayer::import(gtf)
gtf <- as.data.frame(gtf)
refgtf <- as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC148213/sl693/reference/annotation/EMX1_genocde_v40.gtf"))
gtfmerged <- rbind(gtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                   refgtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

counts <- paste0(dir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic_counts.txt")
counts <- fread(counts)
counts <- as.data.frame(counts)

phenotype <- read.csv("/lustre/projects/Research_Project-MRC148213/sl693/SFARI/sfariPhenotype.csv")


## --------------------------- female vs male 

femaleWhole <- phenotype[phenotype$Sex == "F" & grepl("Whole",phenotype$Sample),][["Sample"]]
maleWhole <- phenotype[phenotype$Sex == "M" & grepl("Whole",phenotype$Sample),][["Sample"]]

femaleReads <- counts %>% filter(associated_gene == "EMX1") %>% select(all_of(femaleWhole)) %>% apply(., 1, sum)
maleReads <- counts %>% filter(associated_gene == "EMX1") %>% select(all_of(maleWhole)) %>% apply(., 1, sum)

EMX1_gene <- counts %>% filter(associated_gene == "EMX1") %>% mutate(FReads = femaleReads, MReads = maleReads) %>% filter(exons > 1)


Emx1Iso <- data.frame(
  Isoform = unlist(Emx1Iso <- list(
    Ref = unique(refgtf[refgtf$gene_name == "EMX1" & !is.na(refgtf$transcript_id), "transcript_id"]),
    Female = as.character(EMX1_gene[EMX1_gene$FReads > 0 & EMX1_gene$MReads == 0,][["isoform"]]),
    Male =  as.character(EMX1_gene[EMX1_gene$FReads == 0 & EMX1_gene$MReads > 0,][["isoform"]])

  )),
  Category = rep(names(Emx1Iso), lengths(Emx1Iso))
)
Emx1Iso$colour <- c(rep(NA,nrow(Emx1Iso)))

p = ggTranPlots(inputgtf=gtfmerged, classfiles=counts,
                     isoList = c(as.character(Emx1Iso$Isoform)),
                     selfDf = Emx1Iso, gene="Emx1")


## --------------------------- prenatal vs postnatal 

preWhole <- phenotype[phenotype$Group == "Prenatal" & grepl("Whole",phenotype$Sample),][["Sample"]]
postWhole <- phenotype[phenotype$Group == "Postnatal" & grepl("Whole",phenotype$Sample),][["Sample"]]

preReads <- counts %>% filter(associated_gene == "EMX1") %>% select(all_of(preWhole)) %>% apply(., 1, sum)
postReads <- counts %>% filter(associated_gene == "EMX1") %>% select(all_of(postWhole)) %>% apply(., 1, sum)

EMX1_gene <- counts %>% filter(associated_gene == "EMX1") %>% mutate(preReads = preReads, postReads = postReads) %>% filter(exons > 1)


Emx1Iso <- data.frame(
  Isoform = unlist(Emx1Iso <- list(
    Ref = unique(refgtf[refgtf$gene_name == "EMX1" & !is.na(refgtf$transcript_id), "transcript_id"]),
    Prenatal = as.character(EMX1_gene[EMX1_gene$preReads > 0 & EMX1_gene$postReads == 0,][["isoform"]]),
    Postnatal =  as.character(EMX1_gene[EMX1_gene$preReads == 0 & EMX1_gene$postReads > 0,][["isoform"]])
    
  )),
  Category = rep(names(Emx1Iso), lengths(Emx1Iso))
)
Emx1Iso$colour <- c(rep("#0C0C78",6), rep(NA,nrow(Emx1Iso)-6))

p2 = ggTranPlots(inputgtf=gtfmerged, classfiles=counts,
                isoList = c(as.character(Emx1Iso$Isoform)),
                selfDf = Emx1Iso, gene="Emx1")


## --------------------------- prenatal vs postnatal, male vs female

EMX1_gene <- counts %>% filter(associated_gene == "EMX1") %>% mutate(preReads = preReads, postReads = postReads, FReads = femaleReads, MReads = maleReads) %>% filter(exons > 1)
Emx1Iso <- data.frame(
  Isoform = unlist(Emx1Iso <- list(
    Ref = unique(refgtf[refgtf$gene_name == "EMX1" & !is.na(refgtf$transcript_id), "transcript_id"]),
    Prenatal_female = as.character(EMX1_gene[EMX1_gene$preReads > 0 & EMX1_gene$postReads == 0 & EMX1_gene$FReads > 0 & EMX1_gene$MReads == 0,][["isoform"]]),
    Prenatal_male = as.character(EMX1_gene[EMX1_gene$preReads > 0 & EMX1_gene$postReads == 0 & EMX1_gene$FReads == 0 & EMX1_gene$MReads > 0,][["isoform"]]),
    Postnatal_female =  as.character(EMX1_gene[EMX1_gene$preReads == 0 & EMX1_gene$postReads > 0 & EMX1_gene$FReads > 0 & EMX1_gene$MReads == 0,][["isoform"]]),
    Postnatal_male =  as.character(EMX1_gene[EMX1_gene$preReads == 0 & EMX1_gene$postReads > 0 & EMX1_gene$FReads == 0 & EMX1_gene$MReads > 0,][["isoform"]])
  )),
  Category = rep(names(Emx1Iso), lengths(Emx1Iso))
)
Emx1Iso$colour <- c(rep(NA,nrow(Emx1Iso)))

p3 = ggTranPlots(inputgtf=gtfmerged, classfiles=counts,
                 isoList = c(as.character(Emx1Iso$Isoform)),
                 selfDf = Emx1Iso, gene="Emx1")

plot_grid(p,p2)

nrow(counts[counts$structural_category == "ISM",])/nrow(counts)

