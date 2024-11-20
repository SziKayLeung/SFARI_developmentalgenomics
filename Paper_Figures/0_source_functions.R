#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for sourcing functions for isoform developmental paper
##
## ---------------------------------

suppressMessages(library("dplyr"))
suppressMessages(library("viridis"))
suppressMessages(library("cowplot"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrepel"))
suppressMessages(library("forcats"))
suppressMessages(library("ggh4x"))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(VennDiagram))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

## ---------- Packages -----------------

LOGEN = "C:/Users/sl693/Dropbox/Scripts/LOGen/"
LOGEN_ROOT = "C:/Users/sl693/Dropbox/Scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN, "compare_datasets/whole_vs_targeted.R"))
source(paste0(LOGEN, "compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN, "merge_characterise_dataset/run_ggtranscript.R"))
source(paste0(LOGEN, "differential_analysis/plot_usage.R"))
source(paste0(LOGEN, "aesthetics_basics_plots/draw_venn.R"))
sapply(list.files(path = paste0(LOGEN,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)


## ---------- Labels -----------------

label_colour <- function(genotype){
  if(genotype %in% c("Postnatal")){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "WT_2mos"){colour = alpha(wes_palette("Royal1")[2],0.5)}else{
      if(genotype %in% c("Prenatal")){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "TG_2mos"){colour = alpha(wes_palette("Royal1")[1],0.5)}else{
          if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
            if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
              if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
              }}}}}}}
  return(colour)
}

label_group <- function(genotype){
  if(genotype %in% c("Prenatal")){group = "Prenatal"}else{
    if(genotype %in% c("Postnatal")){group = "Postnatal"}}
  return(group)
}

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


## ---------- CPAT -----------------

plot_cpat <- function(cpat,classfiles){
  
  dat <- merge(cpat[,c("seq_ID","Coding_prob")],classfiles[,c("isoform","structural_category")], by.x = "seq_ID", by.y = "isoform")
  p <- ggplot(dat, aes(x = Coding_prob, y = structural_category, fill = structural_category)) + geom_boxplot() + theme_classic() + 
    labs(x = "Coding probabilty (CPAT)", y = "Structural category") +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values = c("#32cbcf","#b5ebed","#f99189","#fcd6d3","#d9d9d9","#f2ad00","#5d1f1f","#00a08a")) + 
    theme(legend.position = "None")  
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

plot_trans_exp_individual <- function(transcript=NULL, classfiles, Norm_transcounts, var, gene=NULL, sqrt=FALSE, colourdots=NULL){
  
  if(is.null(colourdots)){
    transStru <- classfiles[classfiles$isoform == transcript, "structural_category"]
    structural_colours <- rbind(
      data.frame(
        structural_category = c("Ref","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron","coding","non-coding","noORF", "ORF"),
        structural_col = c("black","#00BFC4",alpha("#00BFC4",0.3),"#F8766D",alpha("#F8766D",0.3),"grey1","gray","grey3","grey4","grey5",wes_palette("Darjeeling1")[2],wes_palette("Royal1")[2],"black","gray")
      )
    )
    colourdots = structural_colours[structural_colours$structural_category == as.character(transStru), "structural_col"]
    print(colourdots)
  }
  
  if(!is.null(transcript)){
    print(transcript)
    dat <- Norm_transcounts %>% filter(isoform == transcript) %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")))
    gene <- classfiles[classfiles$isoform == transcript, "associated_gene"]
  }else{
    print(gene)
    dat <- Norm_transcounts %>% filter(associated_gene == gene) %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal")))
  }
  
  dat <- dat %>% mutate(sex = ifelse(sex == "F", "Female", "Male"))
  print(head(dat))
  
  if(var == "both"){
    p <- ggplot(dat, aes(x = group, y = normalised_counts, colour = transcript)) + 
      geom_boxplot(outlier.shape = NA) + scale_y_sqrt() +
      labs(x = "", y = "Normalized counts", title = paste0(gene, ": ", transcript,"")) +
      scale_colour_manual(values = colourdots) + facet_grid(~sex)
    
  }else{
    if(isTRUE(sqrt)){
      p <- ggplot(dat, aes(x = !!rlang::sym(var), y = normalised_counts, colour = transcript)) + 
        geom_boxplot(outlier.shape = NA) + scale_y_sqrt() +
        labs(x = "", y = "Normalized counts", title = paste0(gene, ": ", transcript,"")) +
        scale_colour_manual(values = colourdots)
    }else{
      p <- ggplot(dat, aes(x = !!rlang::sym(var), y = log10(normalised_counts), fill = !!rlang::sym(var))) + geom_boxplot(outlier.shape = NA) +
        labs(x = "", y = "log10 normalized counts", 
             title = paste0(gene," ", transcript,"")) 
    }
  }

  p <- p + geom_jitter(color=colourdots, size=2, alpha=0.9) +  theme_classic() + 
    theme(strip.background=element_rect(colour="white", fill="white")) +
    #scale_fill_manual(values = c(label_colour(group1),label_colour(group2))) + 
    theme(legend.position = "none") #+ facet_grid(~group) 
  
  
  if(var == "sex"){
    p <- p + labs(x = "Sex") 
  }
  
  return(p)
}

plot_trans_exp_lifetime <- function(transcript=NULL,classfiles,Norm_transcount,gene = NULL,sex=FALSE, colpoints=NULL){
  
  if(is.null(colpoints)){
    transStru <- classfiles[classfiles$isoform == transcript, "structural_category"]
    structural_colours <- rbind(
      data.frame(
        structural_category = c("Ref","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron","coding","non-coding","noORF", "ORF"),
        structural_col = c("black","#00BFC4",alpha("#00BFC4",0.3),"#F8766D",alpha("#F8766D",0.3),"grey1","grey","grey3","grey4","grey5",wes_palette("Darjeeling1")[2],wes_palette("Royal1")[2],"black","gray")
      )
    )
    colpoints = structural_colours[structural_colours$structural_category == transStru, structural_col]
    
  }
  
  #library(grid)
  #library(gridExtra)
  
  if(is.null(gene)){
    gene <- classfiles[classfiles$isoform == transcript, "associated_gene"] 
    dat <- Norm_transcount %>% filter(isoform == transcript) 
    inputTitle <- paste0(gene, ": ", transcript,"")
  }else{
    dat <- Norm_transcount %>% filter(associated_gene == gene) 
    inputTitle <- gene
  }
  
  dat <- dat %>% mutate(group = factor(group, levels = c("Prenatal","Postnatal"))) 

  dat <- as.data.frame(dat)
  fetal.ages.scale <- scales::rescale(dat[which(dat$group=='Prenatal'),'age'], to = c(0, 40))
  child.ages.scale <-  scales::rescale(dat[which((dat$group=='Postnatal') & (dat$age <= 40)),'age'], to = c(41, 60))
  adult.ages.scale <- scales::rescale(dat[which((dat$group=='Postnatal') & (dat$age > 40)),'age'], to = c(75, 100))
  
  dat$age.rescale <- rep(NA, nrow(dat))
  dat$age.rescale[which(dat$group=='Prenatal')] <- fetal.ages.scale
  dat$age.rescale[which((dat$group=='Postnatal')& (dat$age <= 40))] <- child.ages.scale
  dat$age.rescale[which((dat$group=='Postnatal')& (dat$age > 40))] <- adult.ages.scale
  
  FetalBreaksSet <- sort(unique(as.numeric(dat[dat$group == "Prenatal" & dat$age %in% c("6","12","20","28"),"age.rescale"])))
  FetalLabelsSet <- sort(unique(as.numeric(dat[dat$group == "Prenatal" & dat$age %in% c("6","12","20"),"age"])))
  AdultBreaksSet <- sort(unique(as.numeric(dat[dat$group == "Postnatal" & dat$age %in% c("25","41","50","66","80"),"age.rescale"])))
  AdultLabelsSet <- sort(unique(as.numeric(dat[dat$group == "Postnatal" & dat$age %in% c("25","41","50","66","80"),"age"])))
  allBreaksSet <- c(FetalBreaksSet, AdultBreaksSet)
  allLabelsSet <- append(append(FetalLabelsSet,"40/0"),AdultLabelsSet)
  
  dat <<- dat
  if(!isFALSE(sex)){
    p <- ggplot(dat, aes(x = age.rescale, y = log10(normalised_counts), colour = sex)) + geom_point(size = 3) + 
      scale_colour_manual(values = c("pink","blue")) 
  }else{
    p <- ggplot(dat, aes(x = age.rescale, y = log10(normalised_counts))) + geom_point(colour = colpoints, size = 3)
  }
  
  
  p <- p +
    theme_classic() +
    geom_smooth(method=lm, aes(color=isoform), formula = y~poly(x,3),colour="black",fill = alpha("gray",0.2)) +
    theme(panel.spacing = unit(0, "cm", data = NULL),legend.position="top") + 
    labs(x = NULL, y = "log10 normalized counts", title = inputTitle) +
    theme(
      strip.placement = "outside",   # format to look like title
      strip.background = element_blank(),
      legend.position = "None"
    ) +
    force_panelsizes(cols = c(0.4, 1)) +
    #scale_x_continuous(breaks = dat$age.rescale[c(1,15,31,32,33,40,20)], labels = append(append(c(dat$age[c(1,15,31)]),"40/0"), c(dat$age[c(33,40,20)]))) +
    scale_x_continuous(breaks = allBreaksSet, labels = allLabelsSet) +
    geom_vline(xintercept=40, linetype="dotted", color = "black") +
    annotate("text", x = 16, y=max(log10(dat$normalised_counts)) + 0.5, label = "Pre-natal") + 
    annotate("text", x = 70, y=max(log10(dat$normalised_counts)) + 0.5, label = "Post-natal") #+
    #scale_colour_manual(values = c(wes_palette("Royal1")[4],wes_palette("Royal2")[5]))
   
   p <- grid.arrange(p, bottom = textGrob("Age (pcw)                                             Age (yrs)", rot = 0, vjust = 0)) 
   

  
  return(p)
  
}

num_disease_focus_DTE <- function(sigResults, geneList, title=NULL){
  cate_cols <- c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),alpha("#F8766D",0.8),alpha("#F8766D",0.3))
  p <- sigResults %>% 
    mutate(structural_category = factor(structural_category, levels = c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog","genic"))) %>%
    filter(associated_gene %in% geneList) %>% group_by(associated_gene, structural_category) %>% tally() %>% 
    ggplot(., aes(x = reorder(associated_gene,-n), y = n, fill = structural_category)) + geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Target genes", y = "Number of differentially expressed transcripts", title = title) +
    theme(legend.position = c(0.8,0.8)) +
    scale_fill_manual(name = "", values = c(cate_cols,alpha("#808080",0.3)), labels = c("FSM","ISM","NIC","NNC","Genic")) 
  
  return(p)
}

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

## ---------- read lenths pre and post-QC -----------------


plot_lengths <- function(postnatal, prenatal, pfpostnatal, pfprenatal){
  p1<-ggplot(rbind(postnatal, prenatal), aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=rbind(pfpostnatal, pfprenatal),alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Combined\ncount')+scale_y_continuous(labels = scales::scientific)
  p2<-ggplot(prenatal, aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=pfprenatal,alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Prenatal\ncount')+scale_y_continuous(labels = scales::scientific)
  p3<-ggplot(postnatal, aes(x=V2))+geom_density(aes(fill="Pre-filter"),alpha=0.5)+geom_density(aes(fill="Post-filter"), data=pfpostnatal,alpha=0.5)+theme_cowplot()+xlab('Read length')+xlim(-500,65000)+ylab('Postnatal\ncount')+scale_y_continuous(labels = scales::scientific)
  
  return(c(p1,p2,p3))
  
}

## ---------- age distribution -----------------

ages <- function(pheno){
  
  p1<-ggplot(pheno[which(pheno$group=='Postnatal'),], aes(x=age))+geom_histogram(binwidth = 20, colour='black', fill='forestgreen')+ggtitle('Postnatal')+xlab('Age (years)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,15,5))+scale_x_continuous(breaks = seq(0,80,20))+ylab('Number of samples')
  p2<-ggplot(pheno[which(pheno$group=='Prenatal'),], aes(x=age))+geom_histogram(binwidth = 10, colour='black', fill='goldenrod')+ggtitle('Prenatal')+xlab('Age (pcw)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,15,5))+ylab('Number of samples')
  p3<-ggplot(pheno[which(pheno$group=='Postnatal'),], aes(x=age, fill=sex))+geom_histogram(binwidth = 20, position = 'dodge', colour='black')+ggtitle('Postnatal')+xlab('Age (years)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,10,5))+scale_x_continuous(breaks = seq(0,80,20))+ylab('Number of samples')
  p4<-ggplot(pheno[which(pheno$group=='Prenatal'),], aes(x=age, fill=sex))+geom_histogram(binwidth = 10, colour='black', position = 'dodge')+ggtitle('Prenatal')+xlab('Age (pcw)')+theme_cowplot()+scale_y_continuous(limits = c(0,15), breaks=seq(0,10,5))+ylab('Number of samples')
  
  return(list(p1,p2,p3,p4))
         
}


pSensitivity <- function(classf){
  
  dat <- classf %>% arrange(-nreads) %>% mutate(cumreads = cumsum(nreads), relative = prop.table(nreads), cumrel = cumsum(relative)) 
  p <- ggplot(dat, aes(x = isoform, y = cumrel, label = paste0(associated_gene,", ",associated_transcript))) + 
    geom_point() + 
    aes(x = fct_reorder(isoform, cumrel)) + 
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "XX") +
    geom_label_repel(data          = dat[1:5,],
                     size          = 4,
                     box.padding   = 1,
                     point.padding = 0.1,
                     force         = 100,
                     segment.size  = 0.2,
                     segment.color = "grey50",
                     direction     = "x") +
    labs(x = "Transcript", y = "Cumulative read proportion") + mytheme
  
  return(p)
  
  #densityfill <- function(x){
  #  if(x <= 10){return("<10")
  #  }else if (10 < x & x <= 100){return("10-100")
  #  }else if (100 < x & x < 200){return("100-200")
  #  }else {return(">200")}
  #}
  #class.files$targ_all$ndensity <- factor(unlist(lapply(class.files$targ_all$nreads, function(x) densityfill(x))),
  #                                        levels = c("<10","10-100","100-200",">200"))
  #ggplot(class.files$targ_all, aes(x = associated_gene, fill = as.factor(ndensity))) + 
  #  geom_bar() + labs(x = "Target genes", y = "Number of isoforms") + 
  #  scale_fill_discrete(name = "Number of total reads") + mytheme +
  #  theme(legend.position="bottom")
  
}


top_results <- function(diff_results, rank = 10){
  diff_results <- diff_results %>% mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
  top_results <- bind_rows(
    diff_results %>%
      filter(Expression == 'Up-regulated') %>%
      arrange(padj) %>%
      head(rank),
    diff_results %>%
      filter(Expression == 'Down-regulated') %>%
      arrange(padj) %>%
      head(rank)
  )
  return(top_results)
}

plot_top_results <- function(diff_results, exp_results, plot_type, rank=10){
  
  topResults <- top_results(diff_results, rank)
  
  outputDat <- list()
  outputNum = rank * 2
  for(i in 1:outputNum){
    iso = topResults[["isoform"]][i]
    outputDat[[i]] <- plot_trans_exp_individual(iso,exp_results,plot_type)
  }
  upRegOutput = plot_grid(plotlist = outputDat[1:10],ncol=2)
  downRegOutput = plot_grid(plotlist = outputDat[11:20],ncol=2)
  
  retList = list(upRegOutput,downRegOutput)
  names(retList) = c("up","down")
  return(retList)
}

plot_volcano <- function(diff_results,stats=FALSE,interaction="notsex",chromosome=NULL){
  
  #https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/#:~:text=A%20volcano%20plot%20is%20a,tools%20like%20EdgeR%20or%20DESeq2.
  
  if(!is.null(chromosome)){
    # remove X and Y chromosome
    diff_results <- diff_results %>% filter(!chrom %in% c("chrY","chrX"))
  }
  
  diff_results <- diff_results %>% mutate(
      Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                             log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                             TRUE ~ "Unchanged")
  )
  
  message("Number of transcripts:", nrow(diff_results[diff_results$padj < 0.05,]))
  message("Number of transcripts upregulated (red):", nrow(diff_results[diff_results$Expression == "Up-regulated",]))
  message("Number of transcripts downregulated (blue):", nrow(diff_results[diff_results$Expression == "Down-regulated",]))
  

  if(interaction!="sex"){

    Top_genes <- as.data.frame(rbind(diff_results %>% filter(Expression == "Up-regulated") %>% arrange(padj) %>% .[1:10,],
                                     diff_results %>% filter(Expression == "Down-regulated") %>% arrange(padj) %>% .[1:10,]))
    label3 <- c("Pre-natal < Post-natal","Unchanged", "Pre-natal > Post-natal")
    label2 <- c("Pre-natal < Post-natal", "Pre-natal > Post-natal")
    
  }else{
    
    if(is.null(chromosome)){
      Top_genes <- as.data.frame(rbind(diff_results %>% filter(Expression == "Up-regulated") %>% arrange(padj) %>% .[1:10,],
                                       diff_results %>% filter(Expression == "Down-regulated") %>% arrange(padj) %>% .[1:10,]))
    }else{
      Top_genes <- as.data.frame(rbind(diff_results %>% filter(Expression == "Up-regulated" & !chrom %in% c("chrX","chrY")) %>% 
                                         arrange(padj) %>% .[1:10,],
                                       diff_results %>% filter(Expression == "Down-regulated" & !chrom %in% c("chrX","chrY")) %>% 
                                         arrange(padj) %>% .[1:10,]))
    }

    label3 <- c("Female < Male","Unchanged", "Female > Male")
    label2 <- c("Female < Male", "Female > Male")
  }

  
  options(ggrepel.max.overlaps = Inf)
  
  if(isFALSE(stats)){

    diff_results <<- diff_results
    p <- ggplot(diff_results, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion
      geom_point(aes(color = Expression), size = 3/5) +
      xlab(expression("log"[2]*"FC")) +
      ylab(expression("-log"[10]*"FDR")) +
      mytheme +
      guides(colour = guide_legend(override.aes = list(size=2,fill=NA))) +
      theme(legend.key = element_rect(colour = NA, fill = NA)) +
      ggrepel::geom_label_repel(data = Top_genes,
                                mapping = aes(log2FoldChange, -log(padj,10), label = associated_gene),
                                size = 4,
                                box.padding = 1.0,    # Adjust this value to increase/decrease box padding
                                point.padding = 1.0) + theme(legend.position = "top") 
    
    if(nrow(diff_results[diff_results$Expression=="Unchanged",]) == 0){
      p <- p + scale_color_manual(values = c("dodgerblue3", "firebrick3"), label = label2) 
    }else{
      p <- p + scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3"), label = label3) 
    }
    
    output <- list(p, Top_genes)
    names(output) <- c("p","top10")
    return(output)
    
  }
  
}


plotIFTargetedbyGene <- function(gene, pathDIU){
  
  #Exp <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/targeted/group/",gene,"_normalised_expression.txt"),header=T)
  Exp <- read.csv(paste0(pathDIU,"/",gene,"_normalised_expression.txt"),header=T)
  Exp <- Exp %>% tidyr::spread(sample,normalised_counts) %>% tibble::column_to_rownames(var = "isoform")
  p <- plotIF(gene=gene,ExpInput=Exp,pheno=phenotype$WholeTargeted,cfiles=class.files$glob_targ_SQ,design="case_control",rank=5,majorIso=NULL)[[2]]
  return(p)
  
}

plotIFWholebyGene <- function(gene, pathDIU,facetTranscriptsFeature,sexFeature){

  iExp <- fread(pathDIU, data.table = F)
  iExp <- iExp %>% tidyr::spread(sample,normalised_counts) %>% tibble::column_to_rownames(var = "isoform") %>% select(contains("Whole"))
  p <- plotIF(gene=gene,ExpInput=iExp,pheno=phenotype,cfiles=class.files$glob_targ_SQ,design="case_control",rank=5,majorIso=NULL,facetTranscripts=facetTranscriptsFeature,sex=sexFeature)
  return(p)
  
}

fourvenndiagrams <- function(set1, set2,set3,set4,name1, name2, name3,name4){
  p <- venn.diagram(x = list(set1,set2, set3, set4), 
                    label_alpha = 0, category.names = c(name1,name2, name3, name4),filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', 
                    fill = c("#B3E2CD", "#FDCDAC","red","blue"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                    print.mode = "raw")
  return(p)
}

