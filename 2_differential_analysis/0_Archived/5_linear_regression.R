## ---------- Script -----------------
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------



## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/"
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(SC_ROOT, "0_source_differential_functions.R"))
source(paste0(SC_ROOT, "sfari_differential.config.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))

# SFARI gene list
SFARI <- read.csv(paste0(SC_ROOT,"SFARI-Gene_genes_01-23-2023release_01-26-2023export.csv"), header = T)

# Phenotype: pre-natal vs post-natal, male vs female
phenotype$ont_multi$group <- as.character(phenotype$ont_multi$group)
phenotype$ont_multi$time[phenotype$ont_multi$time == 0] <- "prenatal"
phenotype$ont_multi$time[phenotype$ont_multi$time == 1] <- "postnatal"
phenotype$ont_multi$group[phenotype$ont_multi$group == "CASE"] <- "M"
phenotype$ont_multi$group[phenotype$ont_multi$group == "CONTROL"] <- "F"

# Phenotype: fetal, male vs female
colnames(phenotype$brain_whole_fetal) <- c("sample","time","group")

input <- list()
input$phenotype <- phenotype$ont_multi
input$expression <- read.table(talon_expression$brain_whole, sep = "\t", header = T)
input$fetalphenotype <- phenotype$brain_whole_fetal
#input$fetalphenotype_subset <- input$fetalphenotype %>% filter(time < 23)
input$fetalphenotype_subset <- read.table("phenotype.txt", header = T)

# merge expression with classification file
# note only keeping common isoforms in both datasets
#input$expression <- merge(class.files$brain_whole[,c("isoform","associated_transcript","associated_gene")], 
#                          input$expression, by.x = "isoform", by.y = "annot_transcript_id")

# results from running DESeq2
res <- list()
res$prevspostnatal <- run_DESeq2(input$expression,input$phenotype,"annot_transcript_id","M",interaction="On")
res$fetal <- run_DESeq2(input$expression,input$fetalphenotype,"annot_transcript_id","M",interaction="Off")
res$fetalsex <- run_DESeq2(input$fetal_expression,input$fetalphenotype_subset,"annot_transcript_id","M",interaction="On")

# annotate with classification file
col <- c("isoform","associated_transcript","associated_gene")
anno_res <- list()
anno_res$prepost <- merge(res$prevspostnatal$res_Wald, class.files$brain_whole[col], by = "isoform") %>% arrange(padj)
anno_res$fetal <- merge(res$fetal$res_Wald, class.files$brain_whole[,c("isoform","associated_transcript","associated_gene")], by = "isoform") %>% arrange(padj)
anno_res$fetalsex <- merge(res$fetalsex$res_Wald, class.files$brain_whole[,c("isoform","associated_transcript","associated_gene")], by = "isoform") %>% arrange(padj)
anno_res$fetalsexsig <- anno_res$fetalsex %>% filter(anno_res$fetalsex$padj < 0.05, associated_gene %in% SFARI$gene.symbol)
anno_res$fetalsig <- anno_res$fetal %>% filter(anno_res$fetal$padj < 0.05, associated_gene %in% SFARI$gene.symbol)

# plot

p_scatter <- function(transcript, exp_matrix, phenotype, class.file, type){
  print(transcript)
  
  df <- exp_matrix %>%  filter(isoform == transcript) %>% right_join(., phenotype, by = "sample") 
  gene <- class.file[class.file$isoform == transcript,"associated_gene"]


  if(type == "prepost"){

    df$time <- factor(df$time, levels = c("prenatal", "postnatal"))

    p <- ggplot(df, aes(x = time, y = log10(normalised_counts), colour = group)) + 
      geom_boxplot(position=position_dodge2(preserve="single")) + 
      geom_point(aes(fill = group), size = 2.5, shape = 21, position = position_jitterdodge()) +
      labs(x = " ", y = "Normalised counts (Log10)", title = gene, subtitle = transcript) +
      scale_fill_discrete(name = "", labels = c(label_group("Control_timeanalysis"),label_group("Case_timeanalysis"))) +
      scale_x_discrete(drop=FALSE) +
      mytheme +
      guides(color = FALSE, size = FALSE) +
      theme(legend.position = "bottom") 
    
  }else if (type == "fetal"){
    
    df$time <- factor(df$time)
    p <- ggplot(df, aes(x = time, y = log10(normalised_counts), colour = group)) + geom_point() + 
      labs(x = "Age (pcw)", y = "Normalised counts (Log10)", title = gene, subtitle = transcript) +
      scale_colour_discrete(name = "", labels = c(label_group("Control_timeanalysis"),label_group("Case_timeanalysis"))) +
      mytheme +
      theme(legend.position = "bottom") +
      stat_summary(data=df, aes(x=time, y=log10(normalised_counts), group=group), fun ="mean", geom="line", linetype = "dotted", linewidth = 0.5)
    
  }else{
    #df$time <- factor(df$time)
    df <- df %>% filter(time < 23)
    p <- ggplot(df, aes(x = time, y = log10(normalised_counts), colour = group)) + geom_point(size=3) + 
      geom_smooth(method=lm, se=FALSE,linetype = "dotted", linewidth = 1) +
      labs(x = "Age (pcw)", y = "Normalised counts (Log10)", title = gene, subtitle = transcript) +
      scale_colour_discrete(name = "", labels = c(label_group("Control_timeanalysis"),label_group("Case_timeanalysis"))) +
      mytheme +
      theme(legend.position = "bottom") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) 
    
  }
  

  return(p)
}

plot_volcano <- function(diff_results){
  
  #https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/#:~:text=A%20volcano%20plot%20is%20a,tools%20like%20EdgeR%20or%20DESeq2.
  diff_results <- diff_results %>% mutate(
    Expression = case_when(log2FoldChange >= log(2) & `padj` <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & `padj` <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
  
  cat("Number of transcripts upregulated (red):", nrow(diff_results[diff_results$Expression == "Up-regulated",]))
  cat("Number of transcripts downregulated (blue):", nrow(diff_results[diff_results$Expression == "Down-regulated",]))
  
  top <- 10
  top_genes <- bind_rows(
    diff_results %>% 
      filter(Expression == 'Up-regulated') %>% 
      arrange(`padj`, abs(log2FoldChange)) %>% 
      head(top),
    diff_results %>% 
      filter(Expression == 'Down-regulated') %>% 
      arrange(`padj`,abs(log2FoldChange)) %>% 
      head(top)
  )
  
  p <- ggplot(diff_results, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
    geom_point(aes(color = Expression), size = 2/5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    ggrepel::geom_label_repel(data = top_genes,
                              mapping = aes(log2FoldChange, -log(padj,10), label = associated_gene),
                              size = 2)
  
  output <- list(p, top_genes)
  names(output) <- c("p","top10")
  return(output)
  
}

p_scatter("TALONT002306349",res$prevspostnatal$norm_counts,input$phenotype,class.files$brain_whole,"prepost")
p_scatter("ENST00000307845.8",res$fetalsex$norm_counts,input$fetalphenotype[input$fetalphenotype$time < 23,], class.files$brain_whole,"fetal")

ggplot(df, aes(x = time, y = normalised_counts, colour = group)) + geom_point() + geom_smooth(method=lm, se=FALSE)

pPrePost_sex <- lapply(anno_res$prepost$isoform[1:10],function(x) p_scatter(x,res$prevspostnatal$norm_counts,input$phenotype,class.files$brain_whole,"prepost"))
pfetal_nosex <- lapply(anno_res$fetal$isoform[1:10],function(x) p_scatter(x,res$fetal$norm_counts,input$fetalphenotype,class.files$brain_whole,"fetal_nosex"))
pfetal_sex <- lapply(anno_res$fetalsex$isoform[1:10],function(x) p_scatter(x,res$fetalsex$norm_counts,input$fetalphenotype[input$fetalphenotype$time < 23,],class.files$brain_whole,"fetal"))
pfetal_sex_sfari <- lapply(anno_res$fetalsexsig$isoform,function(x) p_scatter(x,res$fetalsex$norm_counts,input$fetalphenotype[input$fetalphenotype$time < 23,],class.files$brain_whole,"fetal"))

pvol <- plot_volcano(anno_res$fetal)
pdf("PrevsPostNatal_top10_DE.pdf", width = 7, height = 5)
for(i in pPrePost_sex){print(i)}
dev.off()

pdf("Fetal_nosexdiff_top10_DE.pdf", width = 7, height = 5)
for(i in pfetal_nosex){print(i)}
dev.off()

pdf("Fetal_sexdiff_top10_DE.pdf", width = 7, height = 5)
for(i in pfetal_sex){print(i)}
dev.off()

pdf("Fetal_sexdiff_sfari_DE.pdf", width = 7, height = 5)
for(i in pfetal_sex_sfari ){print(i)}
dev.off()

pdf("Fetal_nosexdiff_finalposter.pdf", width = 7, height = 5)
p_scatter("TALONT000733011",res$fetal$norm_counts,input$fetalphenotype, class.files$brain_whole,"fetal_nosex")
p_scatter("ENST00000360079.8",res$fetal$norm_counts,input$fetalphenotype, class.files$brain_whole,"fetal_nosex")
dev.off()
res$fetal$res_Wald[res$fetal$res_Wald$isoform == "ENST00000360079.8",]


### manually normalise 
fetal_sample <- intersect(colnames(input$expression), input$fetalphenotype$sample)
input$fetal_expression <- input$expression %>% select("annot_transcript_id", fetal_sample) 
rownames(input$fetal_expression) <- input$fetal_expression$annot_transcript_id
fetal_norm <- apply(input$fetal_expression[-1],2, function(x) x*1000000/sum(x))

fetal_norm_df <- fetal_norm %>% reshape2::melt()
colnames(fetal_norm_df) <- c("isoform","sample","normalised_count")
fetal_norm_df <- merge(fetal_norm_df, input$fetalphenotype, by = "sample")

res$fetalsex <- lm(normalised_count ~ time + group + time*group, data = fetal_norm_df)
#fetal_norm <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/scripts/Normalised_prenatal.csv")
#fetal_norm_df <- reshape2::melt(fetal_norm)

### edgeR and limma for fetal dataset, age vs sex
# https://monashbioinformaticsplatform.github.io/r-linear-abacbs2018/topics/linear_models_abacbs2018.html#5_testing_many_genes_with_limma
library(edgeR)
library(limma)
row.names(input$fetal_expression) <- input$fetal_expression$annot_transcript_id
dgelist <- calcNormFactors(DGEList(input$fetal_expression %>% select(-annot_transcript_id)))
log2_cpms <- cpm(dgelist, log=TRUE, prior.count=0.25)
nrow(log2_cpms)
keep <- rowMeans(log2_cpms) >= -3
log2_cpms_filtered <- log2_cpms[keep,]
nrow(log2_cpms_filtered)

input$fetalphenotype <- input$fetalphenotype[input$fetalphenotype$sample %in% fetal_sample,]
X <- model.matrix(~ time + group + time*group, data=input$fetalphenotype)
fit <- lmFit(log2_cpms_filtered, X)
class(fit)
fit <- eBayes(fit)
topTable(fit, adjust="fdr")
all_results <- topTable(efit, n=Inf)
significant <- all_results$adj.P.Val <= 0.05
table(significant)
ggplot(all_results, aes(x=AveExpr, y=logFC)) + 
  geom_point(size=0.1, color="grey") +
  geom_point(data=all_results[significant,], size=0.1)
