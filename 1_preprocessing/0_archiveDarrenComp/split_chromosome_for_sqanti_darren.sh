#!/bin/bash

# paths
UTILS_DIR=/home/darren/installs/SFARI_developmentalgenomics/1_preprocessing
MERGED_CHROM_DIR=/media/disk2/sfari_RB/5_cupcake/mergedChrom
SPLIT=/media/disk2/sfari_RB/6_sqanti/split

for i in ${MERGED_CHROM_DIR}/*read_stat.txt; do 
  chromosome=$(basename $i .read_stat.txt)
  echo $chromosome
  replaceONTprefix="${chromosome//chr/ONT}"
  #echo $replaceONTprefix
  
  cd ${SPLIT}
  
  # extract the isoform list from each read stat file
  awk -F'\t' '{print $2}' $i | tail -n +2 > WholeTargeted_${chromosome}.id.txt
  
  # replace PB with ONT and the chromosome 
  sed -i "s/PB/${replaceONTprefix}/g" WholeTargeted_${chromosome}.id.txt
  
  # replace gff with the correct gff name
  sed "s/PB/${replaceONTprefix}/g" $MERGED_CHROM_DIR/$chromosome.gff > $MERGED_CHROM_DIR/$chromosome.renamed.gff
  
  # replace read.stat txt with corrected isoform id
  sed "s/PB/${replaceONTprefix}/g" $MERGED_CHROM_DIR/$chromosome.read_stat.txt > $MERGED_CHROM_DIR/$chromosome.read_stat.renamed.txt
  
  # split chunk with 1000000 lines
  split WholeTargeted_${chromosome}.id.txt WholeTargeted_${chromosome}_chunk -l1000000 --additional-suffix=.txt -d

done

ls *chunk* > WholeTargeted_chunks.txt
