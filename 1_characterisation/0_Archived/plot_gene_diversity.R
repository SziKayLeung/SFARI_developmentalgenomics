## ---------- Script -----------------
##
## Purpose of script: Generate plots for SFARI developmental whole transcriptome data (Rosie) 
##
## ---------- Notes -----------------
##
## 
##   
##
##

## ---------- Packages -----------------

library("dplyr")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/sqanti_general.R")


## ---------- Input classification file -----------------

class.files.name = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/8_SQANTI3_RB/sqanti3Filtered2_classification.filtered_lite_classification.txt"
class.files <- SQANTI_class_preparation(class.files.name,"nstandard")

# tally the number of isoforms per associated gene 
transcript_num <- class.files %>% group_by(annot_gene_name) %>% tally()


## ---------- Plots -----------------

# plot the number of isoforms per associated gene by structural category
subset_num_by_gene <- function(gene_list){
  
  p <- class.files %>% filter(associated_gene %in% gene_list) %>%
    group_by(structural_category, associated_gene) %>% tally() %>% 
    ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + 
    labs(x = "Genes", y = "Number of Isoforms") + theme_classic() + 
    theme(legend.position = "bottom") + 
    guides(fill=guide_legend(title="Sturctural Category"))
  
  return(p)
}


# the list of genes to plot
glist = c("MAPT","APP","MEG3")
subset_num_by_gene(glist)