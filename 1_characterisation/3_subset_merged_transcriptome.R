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

