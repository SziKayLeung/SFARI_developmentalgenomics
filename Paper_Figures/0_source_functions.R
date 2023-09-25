#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for sourcing functions for isoform developmental paper
##
## ---------------------------------

suppressMessages(library("viridis"))
suppressMessages(library("cowplot"))
suppressMessages(library("dplyr"))

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


## ---------- read lenths pre and post-QC -----------------


plot_lengths <- function(postnatal, prenatal, pfpostnatal, pfprenatal){
p1<-bind_rows(postnatal, prenatal)%>%filter(V2<10000)%>%ggplot(aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=rbind(pfpostnatal, pfprenatal),alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Combined\ncount')+scale_y_continuous(labels = scales::scientific)
p2<-prenatal%>%filter(V2<10000)%>%ggplot(aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=pfprenatal,alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Prenatal\ncount')+scale_y_continuous(labels = scales::scientific)
p3<-postnatal%>%filter(V2<10000)%>%ggplot(aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=pfpostnatal,alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Postnatal\ncount')+scale_y_continuous(labels = scales::scientific)

return(c(p1,p2,p3))

}

## ---------- age distribution -----------------

ages <- function(pheno){

p1<-ggplot(pheno[which(pheno$group=='Postnatal'),], aes(x=age))+geom_histogram(binwidth = 20, colour='black', fill='forestgreen')+ggtitle('Postnatal')+xlab('Age (years)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,15,5))+scale_x_continuous(breaks = seq(0,80,20))
p2<-ggplot(pheno[which(pheno$group=='Prenatal'),], aes(x=age))+geom_histogram(binwidth = 10, colour='black', fill='goldenrod')+ggtitle('Prenatal')+xlab('Age (pcw)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,15,5))
p2+p1

p3<-ggplot(pheno[which(pheno$group=='Postnatal'),], aes(x=age, fill=sex))+geom_histogram(binwidth = 20, position = 'dodge', colour='black')+ggtitle('Postnatal')+xlab('Age (years)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,10,5))+scale_x_continuous(breaks = seq(0,80,20))
p4<-ggplot(pheno[which(pheno$group=='Prenatal'),], aes(x=age, fill=sex))+geom_histogram(binwidth = 10, colour='black', position = 'dodge')+ggtitle('Prenatal')+xlab('Age (pcw)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,10,5))
p4+p3

return(c(p1,p2,p3,p4)

}

