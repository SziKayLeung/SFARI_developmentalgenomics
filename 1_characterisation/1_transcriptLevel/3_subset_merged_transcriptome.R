library("data.table")
library("dplyr")
demux <-fread("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/demux_fl_count.csv",data.table=FALSE)

demux <- demux %>% filter(isoform != 0)
rownames(demux) <- demux$isoform
demux <- demux %>% select(-merged.fa, -isoform)

demux <- demux %>% mutate(nsamples = rowSums(.!=0), nreads = rowSums(.)) %>% 
  select(nsamples, nreads)

nrow(demux)
oneRead <-  demux %>% filter(nreads == 1 & nsamples == 1)
filtered <- demux %>% filter(nreads >= 2 & nsamples >= 2)
nrow(filtered)

write.table(row.names(filtered),"/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/filteredIds.csv",quote=F, row.names = F, col.names=F, sep = ",")

############## remove intergenic reads
root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
root_sfari <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/"
root_rb_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/"
dirnames <- list(
  # global transcriptome (Iso-Seq, Iso-Seq + RNA-Seq)
  glob_SQ = paste0(root_dir, "/RBFetal/WholeTranscriptome/"),
  
  # targeted sequencing (Iso-Seq, ONT)
  targ_SQ = paste0(root_dir,"RBFetal/1_SQANTI3/"),
  
  wholetarg_SQ = paste0(root_rb_dir,"SQANTI/"),
  
  output = paste0(root_sfari,"0_output/")
)

class.names.files <- list(
  
  # targeted + whole SQANTI dataset futher filtered by minimum 2 reads and 2 samples
  glob_targ_SQ = paste0(dirnames$wholetarg_SQ,"WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification.txt"),
  
  # whole transcriptome SQANTI dataset further filtered by minimum 2 reads and 2 samples
  glob_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples.txt"),
  
  # targeted SQANTI dataset further filtered by minimum 2 reads and 2 samples
  targ_SQ = paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples.txt")
  
) 

class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

# remove mono-exonic intergenic transcripts
class.files <- lapply(class.files, function(x) x %>% mutate(structural_category_exons = paste0(structural_category,"_",exons)))
#mono.class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons == "Intergenic_1"))
class.files <- lapply(class.files, function(x) x %>% filter(structural_category_exons != "Intergenic_1"))

# write the list of isoforms for downstream subsetting
write.table(class.files$glob_targ_SQ$isoform,paste0(dirnames$wholetarg_SQ,"WholeTargeted_RulesFilter_2reads2samples_nomonointergenic_isoforms.txt"), quote=F,row.names=F,col.names = F)
write.table(class.files$glob_SQ$isoform,paste0(dirnames$wholetarg_SQ,"WholeTargeted_RulesFilter_Whole_2reads2samples_nomonointergenic_isoforms.txt"), quote=F,row.names=F,col.names = F)

