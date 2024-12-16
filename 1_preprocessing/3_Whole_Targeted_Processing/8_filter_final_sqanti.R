#!/usr/bin/env Rscript
## ----------Script-----------------
## 
## SK.Leung: further filter by 2 reads, 2 samples
## merge sqanti filtered by mono-exonic transcripts with expression file (demux)
## ---------------------------------

library("data.table")
library("dplyr")

dirnames <- list(
  	output ="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/",
  	readstat = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux/2_finalised_readstat",
  	demux = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux/4_demux_merged"
)


# read in sqanti classification file (whole + targeted merged across 87 datasets)
# sqanti filtered (relaxed json file), removal of mono-exonic intergenic transcripts, removal of mono-exonic transcripts within multi-exonic genes
class.files <- fread(paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_classification.txt"), sep = "\t", data.table=FALSE)

# read in demux files across 24 chromosomes  
demux.names.files <- list.files(path = dirnames$demux, pattern = "fl_count.csv", full = TRUE)
demux.files <- lapply(demux.names.files, function(x) fread(x, data.table=FALSE))
merged.demux.files <- bind_rows(demux.files)

# read stat with minimum 2 FL reads 
readstat.names.files <- list.files(path = dirnames$readstat, pattern = "min2FL.sqantiMonoexonicFilteredID.txt", full = TRUE)
readstat.files <- lapply(readstat.names.files, function(x) fread(x, data.table=FALSE, header=FALSE))
readstatID <- unlist(lapply(readstat.files, function(df) df[[1]]))

sqantikept_2FL <- intersect(as.character(readstatID),as.character(class.files$isoform))

if(length(setdiff(sqantikept_2FL, as.character(merged.demux.files$id))) != 0){print("Error: mismatched IDs")}
if(length(setdiff(as.character(merged.demux.files$id), as.character(readstatID))) != 0){print("Error: missing demux IDs")}
# "pbid" in readstatID as extracted the second column of read.stat.files which includes header, but not in merged.demux.files 
if(length(setdiff(as.character(readstatID), as.character(merged.demux.files$id))) > 1){print("Error: missing demux IDs")}

# check all merged.demux.files nreads >= 2
# demux generated from read_stat_min2FL
if(nrow(merged.demux.files %>% filter(nreads < 2)) != 0){print("Error: other transcripts retained with less than 2 read count support")}

# filter demux files
filtered.merged.demux.files <- merged.demux.files %>% filter(nreads >= 2 & nsamples >= 2)

filtered.class.files <- class.files %>% filter(isoform %in% filtered.merged.demux.files$id)

# final classification file 
final.class.files <- merge(filtered.class.files, filtered.merged.demux.files, by.x = "isoform", by.y = "id", all.x = TRUE)

write.table(final.class.files, file = paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt"), quote=F, sep = "\t", row.names= F)
write.table(final.class.files$isoform, file = paste0(dirnames$output,"sqantifiltered_monoexonic_2reads2samplesfiltered_ID.txt"), quote=F, sep = "\t", row.names= F, col.names = F)

### ------ separate filter by expression: 10 reads, 10 samples ----- 

filtered.merged.demux.files <- merged.demux.files %>% filter(nreads >= 10 & nsamples >= 10)
filtered.class.files <- class.files %>% filter(isoform %in% filtered.merged.demux.files$id)

# final classification file 
final.class.files <- merge(filtered.class.files, filtered.merged.demux.files, by.x = "isoform", by.y = "id", all.x = TRUE)

write.table(final.class.files, file = paste0(dirnames$output,"sqantifiltered_monoexonicfiltered_10reads10samples_classification.txt"), quote=F, sep = "\t", row.names= F)
write.table(final.class.files$isoform, file = paste0(dirnames$output,"sqantifiltered_monoexonic_10reads10samplesfiltered_ID.txt"), quote=F, sep = "\t", row.names= F, col.names = F)