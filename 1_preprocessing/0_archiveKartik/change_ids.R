library(tidyverse)
library(data.table)

for(i in seq(1, 22)){
x<-fread(paste0('Whole_cleaned_aligned_merged_collapsed_chr',i,'.read_stat.txt'), stringsAsFactors=F, data.table=F)

x$pbid<-str_replace_all(string=str_replace(string=x$pbid, pattern='PB', replacement=paste0('ONT', i)), pattern='\\.', replacement='_')

fwrite(x, paste0('Whole_cleaned_aligned_merged_collapsed_chr',i,'.read_stat.renamed.txt'), row.names=F, quote=F, sep='\t')

}
