library("dplyr")
library("stringr")

## -------------- directory paths and input files ----------------

output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/0_flask/1_metadata/"

metadata = read.table("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/10_characterisation_SK/0_metadata//Brain_TargetedTranscriptome_SFARI_metadata.txt", header = T) 

demux = data.table::fread("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/FINAL/demux_fl_count.csv", sep = ",")


## ------------- discard isoforms < 5 reads --------------

# remove .fa from column names
colnames(demux) <- word(colnames(demux), c(1), sep = fixed("."))

# remove isoform == 0, left over isoforms that were not collapsed but discarded
demux <- demux %>% filter(isoform != "0")

# sum the number of reads across each isoform and filter isoform < 5 reads
demuxSum = demux[,-c("isoform")] %>% apply(., 1, sum)
demux$FLRead = demuxSum
demuxFiltered <- demux %>% filter(demux$FLRead >= 5)


## ------------- output --------------

# demux file 
demuxFiltered <- demuxFiltered %>% select("isoform", intersect(colnames(demux), metadata$SQ_ID))
write.csv(demuxFiltered, paste0(output_dir, "demux.csv"), quote = F, row.names = F)

# metadata file
metadata <- metadata %>% filter(SQ_ID %in% intersect(colnames(demux), metadata$SQ_ID))
write.csv(metadata, paste0(output_dir, "phenotype.csv"), quote = F, row.names = F)
