#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address

# 10/02/2024: subset bam file for HNRNPK reads for IGV 

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
grep_raw_reads ONT9.2747.18309 chr9
grep_raw_reads ONT9.2747.18433 chr9
grep_raw_reads ONT9.2747.20508 chr9
grep_raw_reads ONT9.2747.10079 chr9
grep_raw_reads ONT9.2747.18418 chr9
grep_raw_reads ONT1.1323.1553 chr1
grep_raw_reads ONT10.834.859 chr10
grep_raw_reads ONT21.2978.4967 chr21
