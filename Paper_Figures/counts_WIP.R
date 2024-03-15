#counts = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/Targeted_FINAL_demux_2reads2samples.csv")
counts = data.table::fread("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_demux_2reads2samples_SQANTIfiltered.csv")

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/compare_datasets/dataset_identifer.R"))


dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/",
  root_sfari = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/",
  wholetarg_SQ = "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/",
  testing = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/"
)
phenotype <- list(
  WholeTargeted = read.csv(paste0(dirnames$root_sfari, "/12_deseq2/WholeTargetedphenotype.csv")) %>% mutate(time = age)
)
phenotype$WholeTargeted <- phenotype$WholeTargeted %>% mutate(col = paste0(sample,"_",group))

# replace column names of rawExp with mathing col (with group) 
names(counts) <- phenotype$WholeTargeted$col[match(names(counts), phenotype$WholeTargeted$sample)]
names(counts)[1] <- "isoform"
names(counts)[length(counts)-1] <- "nreads"
names(counts)[length(counts)] <- "nsamples"

counts$prenatal_sum_FL <- apply(counts %>% select(contains("Prenatal")),1,sum)
counts$postnatal_sum_FL <- apply(counts %>% select(contains("Postnatal")),1,sum)

counts$dataset <- apply(counts, 1, function(x) identify_dataset_by_counts(x[["prenatal_sum_FL"]], x[["postnatal_sum_FL"]], "prenatal","postnatal"))
write.csv(counts, "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_demux_2reads2samples_SQANTIfiltered_dataset.csv",
          quote=F,row.names=F)
