#!/bin/bash

# 10/02/2024: subset bam file for HNRNPK reads for IGV 

module load SAMtools
Collapse=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
Junctions=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_junctions.txt

cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/IGV

grep_raw_reads(){
	
	transcript=$1
	echo "Processing $transcript"
	grep ${transcript} ${Collapse}/chr9.read_stat.renamed.min2FL.txt | awk '{print $1}' > ${transcript}.RawReads.txt
	samtools view -h ${Collapse}/chr9.bam | grep -f ${transcript}.RawReads.txt > body.sam
	cat chr9header.sam body.sam > ${transcript}.RawReads.sam
	echo "Sam file generated"
	rm body.sam
	samtools view -bS ${transcript}.RawReads.sam > ${transcript}.RawReads.bam
	samtools index ${transcript}.RawReads.bam
	
	grep $transcript ${Junctions}
	
}
grep_raw_reads ONT9.2747.18309
grep_raw_reads ONT9.2747.18433
grep_raw_reads ONT9.2747.20508	
grep_raw_reads ONT9.2747.10079
grep_raw_reads ONT9.2747.18418

