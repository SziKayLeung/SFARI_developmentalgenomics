#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for sourcing functions for isoform developmental paper
##
## ---------------------------------

## ---------- Packages -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))


## ---------- numIso -----------------

# Aim: plot the number of isoforms
# Pre-requisite: SQANTI_gene_preparation() from LOGEN/read_sq_classification.R
# Input:
  # classification file generated from SQANTI
# Output:
  # p = bar plot of the percentage of isoforms per gene

numIso <- function(class.files){
  
  p <- SQANTI_gene_preparation(class.files) %>%
    mutate(cate = ifelse(structural_category %in% c("FSM","ISM"),"Known","Novel")) %>%
    group_by(nIsoCat, cate) %>% tally(nIso) %>%
    mutate(Perc = n/sum(n) * 100) %>%
    ggplot(., aes(x=nIsoCat,y= Perc, fill = cate)) +
    geom_bar(stat="identity", position = position_dodge()) + 
    labs(x ="Number of Isoforms", y = "Genes (%)", fill = "", title = "\n") +
    theme_classic()
  
  return(p)
}


## ---------- numIsoCate -----------------

# Aim: plot the number of isoforms by structural category
# Pre-requisite: tabulate_structural_cate() from LOGEN/read_sq_classification.R
# Input:
  # classification file generated from SQANTI
# Output:
  # p = bar plot of the percentage of isoforms by structural category

numIsoCate <- function(class.files){
  
  p <- tabulate_structural_cate(class.files) %>%
    filter(!is.na(structural_category)) %>%
    ggplot(., aes(x = structural_category, y = perc)) + 
    geom_bar(stat="identity") + labs(x = "Structural Category", y = "Isoforms (%)") +
    coord_flip() 
  
  return(p)
  
}


## ---------- targetRate -----------------

targetRate <- function(){
  # on-target rate
  expressionAll %>% rownames_to_column(., var = "isoform") %>% 
    mutate(target = ifelse(isoform %in% ClassFiles$targ$isoform,"Target","OffTarget")) %>% 
    filter(isoform != "0") %>%
    group_by(target) %>% tally(total) %>% mutate(perc = n/sum(n) * 100)
  
  
  # off target expression
  offTargetExp <- expressionAll %>% rownames_to_column(., var = "isoform") %>% mutate(target = ifelse(isoform %in% ClassFiles$targ$isoform,"Target","OffTarget")) %>% 
    filter(target == "OffTarget") %>%
    left_join(., ClassFiles$all[,c("isoform","associated_gene")], by = "isoform") %>%
    filter(isoform != "0") %>%
    group_by(associated_gene) %>% tally(total)
  
  View(offTargetExp)
}
