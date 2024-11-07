cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/13_massSpec
find . -type f -name *_novel_peptides.gtf* -exec cat {} + > AllSfari_novelpeptides.gtf
find . -type f -name '*novel_peptides_seq.bed12*' -exec cat {} + > Allsfari_novel_peptides.bed12
sort Allsfari_novel_peptides.bed12 | uniq > Allsfari_novel_peptides_unique.bed12
sort AllSfari_novelpeptides.gtf | uniq > AllSfari_novelpeptides_unique.gtf