#!/usr/bin/env Rscript
## ----------Script-----------------
## 
## SK.Leung: generate list of ID of isoforms detected in the whole transcriptome dataset 
## minimum 2 reads, 2 samples
## not detected in the targeted dataset
## using K.Chundruv final version classification text modified for novelGene name nomenclature only 
## ---------------------------------

library("data.table")
library("dplyr")

dirnames <- list(
	sqanti = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/"
)

# read in classification file
class.files <- fread(paste0(dirnames$sqanti, "sqantifiltered_monoexonicfiltered_2reads2samples_classification_finalversion.txt"), data.table = FALSE)

# filter isoforms that are only detected in the whole dataset, 2 reads, 2 samples
whole.class.files <- class.files %>% 
	filter(targeted_nsamples == 0 , targeted_nreads == 0) %>% 
	filter(whole_nsamples >= 2, whole_nreads >= 2)
write.table(whole.class.files$isoform, paste0(dirnames$sqanti,"sqantifiltered_monoexonic_2reads2samplesfiltered_whole_ID.txt"), quote=F, sep = "\t", row.names= F, col.names = F)

