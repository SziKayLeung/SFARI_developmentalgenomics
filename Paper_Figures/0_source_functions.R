#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for sourcing functions for isoform developmental paper
##
## ---------------------------------

suppressMessages(library("viridis"))
suppressMessages(library("cowplot"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrepel"))
suppressMessages(library("forcats"))

## ---------- Packages -----------------

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "compare_datasets/whole_vs_targeted.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))


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

plot_trans_exp_individual <- function(transcript, Norm_transcounts, var){
  print(transcript)
  dat <- Norm_transcounts %>% filter(isoform == transcript)
  gene <- dat$associated_gene[1]
  
  p <- ggplot(dat, aes(x = !!rlang::sym(var), y = log10(normalised_counts), fill = !!rlang::sym(var))) + geom_boxplot(outlier.shape = NA) + 
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    labs(x = "", y = "Isoform Expression (log10)", 
         title = paste0(gene, ":", transcript,"")) +  theme_classic() + 
    #scale_fill_manual(values = c(label_colour(group1),label_colour(group2))) + 
    theme(legend.position = "none") #+ facet_grid(~group)
  
  return(p)
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

plot_volcano <- function(diff_results,stats=FALSE){
  
  #https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/#:~:text=A%20volcano%20plot%20is%20a,tools%20like%20EdgeR%20or%20DESeq2.
  diff_results <- diff_results %>% mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
  
  message("Number of transcripts:", nrow(diff_results[diff_results$padj < 0.05,]))
  message("Number of transcripts upregulated (red):", nrow(diff_results[diff_results$Expression == "Up-regulated",]))
  message("Number of transcripts downregulated (blue):", nrow(diff_results[diff_results$Expression == "Down-regulated",]))
  
  if(isFALSE(stats)){

    
    p <- ggplot(diff_results, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion
      geom_point(aes(color = Expression), size = 2/5) +
      xlab(expression("log"[2]*"FC")) +
      ylab(expression("-log"[10]*"FDR")) +
      scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
      guides(colour = guide_legend(override.aes = list(size=1.5))) + mytheme  +
      ggrepel::geom_label_repel(data = top_genes,
                                mapping = aes(log2FoldChange, -log(padj,10), label = associated_gene),
                                size = 2) + theme(legend.position = "top")
    
    output <- list(p, top_genes)
    names(output) <- c("p","top10")
    return(output)
    
  }
  
}

