#!/bin/bash

# paths
pbmm2Aligned=/media/disk2/sfari_RB/5_cupcake/5_align
TargetedSamples=($(ls ${pbmm2Aligned}/*mapped_filtered_sorted*bam))
WholeSamples=($(ls ${pbmm2Aligned}/*aligned_clean_aligned_filtered.bam))

mamba activate SQANTI3.env

cd ${pbmm2Aligned}
for i in ${TargetedSamples[@]}; do 
	echo "Processing: $i"
	sample=$(basename $i .bam)
	echo $sample
	samtools bam2fq $i | seqtk seq -A > $sample.fa
done

for i in ${WholeSamples[@]}; do 
	echo "Processing: $i"
	sample=$(basename $i .bam)
	echo $sample
	samtools bam2fq $i | seqtk seq -A > $sample.fa
done
