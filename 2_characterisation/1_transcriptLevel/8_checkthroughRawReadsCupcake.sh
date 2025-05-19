#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=8_checkthroughRawReadsCupcake.o
#SBATCH --error=8_checkthroughRawReadsCupcake.e

cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/rawCollapsed 

# 1. from SQANTI classification file, obtain the isoform ID associated with HNRPNK
#sqantiClassFile=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt
#grep HNRNPK ${sqantiClassFile} | awk -F ' ' '{print $2}' > HNRNPK_finalID.txt

# 2. use finalID to obtain the original raw ONT read ID that it was collapsed from
#collapsedFile=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/WholeTargeted_cleaned_aligned_merged_collapsed.read_stat.renamed.txt
#grep -w -f HNRNPK_finalID.txt ${collapsedFile} | awk -F ' ' '{print $1}'  > HNRNPK_rawONTID.txt

# 3. merge all sorted.bam files from aligned minimap2
# note only 39 samples, as missing samples from transfer
module load Miniconda2
source activate nanopore

alignDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/3_aligned
for i in ${alignDir}/*sorted.bam; do
  sample=$(basename "$i" | cut -d "." -f 1 )
  echo $sample
  samtools view $i | grep -w -f HNRNPK_rawONTID.txt > ${sample}_HNRNPK.sorted.sam
  samtools view -bS ${sample}_HNRNPK.sorted.sam > ${sample}_HNRNPK.sorted.bam
  samtools index  ${sample}_HNRNPK.sorted.bam
done

samtools merge mergedWholeTargeted_HNRNPK.sorted.sam *HNRNPK.sorted.sam

# 4. extract raw ONT read ID 
samtools merge mergedWholeTargeted.bam ${alignDir}/*sorted.bam
samtools view mergedWholeTargeted.bam | grep -f HNRNPK_rawONTID.txt > mergedWholeTargeted_HNRNPK.sorted.sam
samtools view -bS mergedWholeTargeted_HNRNPK.sorted.sam > mergedWholeTargeted_HNRNPK.sorted.bam
samtools index mergedWholeTargeted_HNRNPK.sorted.bam mergedWholeTargeted_HNRNPK.sorted.bam.bai