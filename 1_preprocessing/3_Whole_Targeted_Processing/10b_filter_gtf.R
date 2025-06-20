library(data.table)
library(plyranges)

y<-fread('/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_classification_finalversion.txt', stringsAsFactors=F, data.table=F)
for(i in c(seq(1,22), 'X', 'Y')){
 x<-read_gff(paste0('scratch/chr',i,'.gff'))
 x<-x[which(x$transcript_id%in%y$isoform),]
 write_gff3(x, paste0('chr',i,'.gff'))}

