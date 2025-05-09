#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Additional plots after 2nd revision (May 2025)
## 	sensitivity curve of number of reads and samples in whole dataset only
## --------------------------------

## ---------- packages -----------------

library("data.table")

## ---------- paths -----------------

LOGenDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGenDir,"transcriptome_stats/sample_sensitivity.R"))
source(paste0(LOGenDir, "aesthetics_basics_plots/pthemes.R"))

sqanti <- fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron_classification_updatednames.txt",
data.table = F)

## ---------- sensitivity analysis -----------------

sqanti$cate_sumFL <- categories(as.numeric(as.character(sqanti$whole_nreads)))
dat2 <- sqanti %>%
  group_by(whole_nsamples, cate_sumFL) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(cate_sumFL = factor(cate_sumFL, levels = c(seq(1:9),"10","10-20","21-30","31-40","41-50","> 50")))
ggplot(dat2, aes(x = as.factor(whole_nsamples), y = proportion, colour = cate_sumFL, group = cate_sumFL)) +
  geom_point() +
  geom_line() + scale_y_continuous(labels = perc_lab)  + 
  labs(x = "Minimum number of samples", y = "Reads (%)",
       colour = "Minimum number of reads") +
  theme_classic() 

# 10 reads, 10 samples
nrow(sqanti %>% filter(whole_nreads >= 10 & whole_nsamples >= 10))

# 20 reads, 12 samples
nrow(sqanti %>% filter(whole_nreads >= 20 & whole_nsamples >= 12))

# type of isoforms filtered after applying a more stringent filter
sqantiFurtherFiltered <- sqanti %>% filter(whole_nreads >= 20 & whole_nsamples >= 12)
diff <- sqanti[!sqanti$isoform %in% sqantiFurtherFiltered$isoform,]
diff %>% group_by(structural_category, subcategory, monomulti) %>% tally() %>%
ggplot(., aes(x = structural_category, y = n, fill = subcategory)) + 
	geom_bar(stat = "identity") + facet_grid(~monomulti) +
	labs(y = "number of transcripts removed")