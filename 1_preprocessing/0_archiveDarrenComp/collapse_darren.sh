## --------------------------------------
#!/bin/bash

# paths
UTILS_DIR=/home/darren/installs/SFARI_developmentalgenomics/1_preprocessing
MERGED_CHROM_DIR=/media/disk2/sfari_RB/5_cupcake/mergedChrom

# list the pbmm2 aligned individual files to be merged
# create file if not present
if [ ! -f ${UTILS_DIR}/combined_files_darren.txt ]; then 
  ls /media/disk2/sfari_RB/5_cupcake/5_align/*Whole*filtered.bam > ${UTILS_DIR}/combined_files_darren.txt
  ls /media/disk2/sfari_RB/5_cupcake/5_align/*Targeted*filtered_sorted.bam >> ${UTILS_DIR}/combined_files_darren.txt
fi

mamba activate SQANTI3.env
##-------------------------------------------------------------------------

# merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y)
cd ${MERGED_CHROM_DIR}
echo "Merging"
samtools merge -f WholeTargeted.bam -b ${UTILS_DIR}/combined_files_darren.txt
echo "sorting"
samtools sort WholeTargeted.bam -o WholeTargeted_sorted.bam 
samtools index ${MERGED_CHROM_DIR}/WholeTargeted_sorted.bam

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
for i in {0..24}; do 
	echo $i	
	chrNum=${chromosomes[${i}]}  
	echo $chrNum
	
	# isoseq collapse the merged file by chromosome
	samtools view -b ${MERGED_CHROM_DIR}/WholeTargeted_sorted.bam $chrNum > ${MERGED_CHROM_DIR}/${chrNum}.bam

	isoseq3 collapse ${MERGED_CHROM_DIR}/${chrNum}.bam ${chrNum}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16

done
