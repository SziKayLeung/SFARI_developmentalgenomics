module load SAMtools
Collapse=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
Junctions=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_junctions.txt

cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/IGV

grep_raw_reads(){
	
	transcript=$1
 	chr=$2
  
	echo "Processing $transcript"
	samtools view -H ${Collapse}/${chr}.bam > ${chr}header.sam
 	grep ${transcript} ${Collapse}/${chr}.read_stat.renamed.min2FL.txt | awk '{print $1}' > ${transcript}.RawReads.txt
	samtools view -h ${Collapse}/${chr}.bam | grep -f ${transcript}.RawReads.txt > ${chr}body.sam
	cat ${chr}header.sam ${chr}body.sam > ${transcript}.RawReads.sam
	echo "Sam file generated"
	#rm ${chr}body.sam
	samtools view -bS ${transcript}.RawReads.sam > ${transcript}.RawReads.bam
	samtools index ${transcript}.RawReads.bam
	
	grep $transcript ${Junctions}
	
}
grep_raw_reads ONT8.85.5521 chr8
grep_raw_reads ONT3.117.12367​ chr3

SQANTIDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/
GTF=${SQANTIDIR}/sqantifiltered_monoexonicfiltered_2reads2samples_intergenicGenicIntron.filtered.gtf
grep -w ONT8.85.5521 ${GTF} > ${SQANTIDIR}/ONT8.85.5521.gtf
grep -w ONT3.117.12367 ${GTF} > ${SQANTIDIR}/ONT3.117.12367.gtf