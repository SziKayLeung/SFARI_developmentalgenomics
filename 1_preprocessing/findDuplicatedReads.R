library("dplyr")
library("ggplot2")
library("cowplot")
library("data.table")

dir <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/2_trimmed/targeted/duplicated/"
targetedTotalNumber <- read.csv(paste0(dir, "targeted_totalNumber_Reads.csv"), header = F, sep = "")
targetedDuplicated <- read.csv(paste0(dir, "targeted_totalNumber_duplicatedReads.csv"), header = F, sep = "")

manifest <- read.csv("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/WholeTargetedphenotype_manifest.csv")

Numbers <- merge(targetedTotalNumber, targetedDuplicated , by = "V2")
colnames(Numbers) <- c("ID","totalNumbers","duplicatedNumbers")
Numbers <- Numbers %>% mutate(perc = duplicatedNumbers/totalNumbers * 100) 
Numbers <- merge(Numbers, manifest, by = "ID")

p1 <- ggplot(Numbers, aes(x = RIN, y = totalNumbers, colour = Group)) + geom_point(size = 3) + theme_classic() + scale_y_continuous(trans='log10') +
  geom_smooth(method='lm', formula= y~x,colour="black")
p2 <- ggplot(Numbers, aes(x = RIN, y = perc)) + geom_point(size = 3) + theme_classic() + 
  labs(y = "Percentage of duplicated reads/total reads", x = "RIN") + 
  geom_smooth(method='lm', formula= y~x,colour="black")
plot_grid(p1,p2, labels = c("A","B"))

cor.test(Numbers$totalNumbers, Numbers$RIN)
cor.test(Numbers$duplicatedNumbers, Numbers$RIN)


## ------ 

duplicatedReadNames <- list.files(path = dir, pattern = "_duplicated_read_stats_final.csv", full.names = T)
duplicatedReads <- lapply(duplicatedReadNames, function(x) fread(x, header = F))
names(duplicatedReads) <- stringr::word(basename(duplicatedReadNames),c(6),sep=stringr::fixed("_"))

duplicatedReadsMerged <- do.call(rbind, duplicatedReads)
names(duplicatedReadsMerged) <- c("ONTRawReadID","CollapsedID")

#AlldemuxReads <- WholeTargeted_demux.csv
load(file = paste0(root_rb_dir, "6_sqanti3/all_filtered_classification_2reads2samples_noMonoIntergenicAll.RData"))

message("Number of isoforms: ", nrow(class.files$glob_targ_SQ))
message("Number of duplicated reads: ", length(duplicatedReadsMerged$CollapsedID))
duplicatedClass.files <- class.files$glob_targ_SQ %>% filter(isoform %in% duplicatedReadsMerged$CollapsedID)
message("Number of duplicated reads in final dataset: ", nrow(duplicatedClass.files))
table(duplicatedClass.files$structural_category)
duplicatedReadsMerged <- as.data.frame(duplicatedReadsMerged)

finalDuplicatedNumber <- duplicatedReadsMerged[duplicatedReadsMerged$CollapsedID %in% duplicatedClass.files$isoform, ] %>% group_by(ONTRawReadID) %>% 
  tally()

duplicatedReadsMerged[duplicatedReadsMerged$ONTRawReadID == "e6e44c22-4720-40ac-9746-46bf0ca0f9c2",]

duplicatedReadsMerged[duplicatedReadsMerged$ONTRawReadID == "00000070-91e1-494f-85be-132a12b8164d",]
duplicatedReadsMerged[duplicatedReadsMerged$CollapsedID == "ONT3_7072_4054",]
duplicatedReadsMerged[duplicatedReadsMerged$CollapsedID == "ONT3_7070_4470",]

duplicatedClass.files[duplicatedClass.files$structural_category == "NIC",]
