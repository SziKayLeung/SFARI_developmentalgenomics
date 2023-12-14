#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: code for descriptive stats for isoform developmental paper
##
## ---------------------------------

# whole transcriptome
numIsogene = numIsoGene(class.files$glob_SQ,stats=TRUE)
numIsogeneTally = class.files$glob_SQ %>% group_by(associated_gene) %>% tally()
min(numIsogeneTally$n)
median(numIsogeneTally$n)
max(numIsogeneTally$n)
sd(numIsogeneTally$n)
nrow(numIsogeneTally[numIsogeneTally$n > 1,])
nrow(numIsogeneTally[numIsogeneTally$n > 1,])/length(unique(numIsogeneTally$associated_gene))

nrow(numIsogeneTally[numIsogeneTally$n >= 10,])
nrow(numIsogeneTally[numIsogeneTally$n >= 10,])/length(unique(numIsogeneTally$associated_gene))

nrow(mono.class.files$glob_SQ)/nrow(class.files$glob_SQ)

message("prenatal vs postnatal")
message("Number of differentially expressed isoforms: ", nrow(WholeDESeqSig$age))
message("Number of genes with differentially expressed isoforms: ", length(unique(WholeDESeqSig$age$associated_gene)))
message("Number of differentially expressed FSM isoforms: ", 
        nrow(WholeDESeqSig$age %>% filter(structural_category %in% c("FSM"))))
message("Number of known genes with differentially expressed isoforms: ", 
        length(unique(WholeDESeqSig$age$associated_gene[!grepl("novel", WholeDESeqSig$age$associated_gene)])))


plot_volcano(diff_results=WholeDESeq2$age,stats = TRUE)
nrow(WholeDESeq2Sig$age %>% group_by(associated_gene) %>% tally())
WholeDESeqSig$age %>% group_by(structural_category) %>% tally()
17308/85428
63204/85428

nrow(WholeDESeq2Sig$sex)
nrow(WholeDESeq2Sig$sex %>% filter(!chrom %in% c("chrX","chrY")))

message("Mean number of exons, (sd): ", round(mean(class.files$glob_SQ$exons),2)," (", round(sd(class.files$glob_SQ$exons),2),")")
class.files$targ_SQ %>% group_by(associated_gene) %>% tally() %>% arrange(-n)

Cpat$whole

# targeted transcriptome
message("Mean number of exons, (sd): ", round(mean(class.files$targ_SQ$exons),2)," (", round(sd(class.files$targ_SQ$exons),2),")")
length(unique(class.files$targ_SQ$associated_gene))
nrow(class.files$targ_SQ)
mean(class.files$targ_SQ$length)
sd(class.files$targ_SQ$length)
min(class.files$targ_SQ$length)
max(class.files$targ_SQ$length)
numIsogeneTallyTargeted = class.files$targ_SQ%>% group_by(associated_gene) %>% tally()
min(numIsogeneTallyTargeted$n)
max(numIsogeneTallyTargeted$n)
median(numIsogeneTallyTargeted$n)
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 1,])
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 1,])/length(unique(numIsogeneTallyTargeted$associated_gene))
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 10,])
nrow(numIsogeneTallyTargeted[numIsogeneTallyTargeted$n > 10,])/length(unique(numIsogeneTallyTargeted$associated_gene))
nrow(class.files$targ_SQ[class.files$targ_SQ$structural_category == "NNC",])/nrow(class.files$targ_SQ)

plot_volcano(diff_results=TargetedDESeq2$age,stats = TRUE)
nrow(TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n>10))
TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% arrange(-n)
159/length(unique(TargetedDESeq2Sig$age$associated_gene))
TargetedDESeq2Sig$age %>% group_by(structural_category) %>% tally()

TargetedDESeq2Sig$age %>% filter(associated_gene %in% c(TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n==1) %>% .$associated_gene))
TargetedDESeq2Sig$age %>% group_by(associated_gene) %>% tally() %>% filter(n==1)
807/nrow(TargetedDESeq2Sig$age)
4692/nrow(TargetedDESeq2Sig$age)
