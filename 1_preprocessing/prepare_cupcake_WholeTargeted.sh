#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=prepare_cupcake_WholeTargeted.o
#SBATCH --error=prepare_cupcake_WholeTargeted.e


echo Job started on:
date -u

module load Miniconda2/4.3.21

# paths
UTILS_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/0_utilities
ALIGNED_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/pbmm2_align
MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
MERGED_CUPCAKE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/cupcakeMerged

# sort whole.bam
source activate nanopore
#for i in ${ALIGNED_WHOLE_DIR}/*Whole*filtered.bam; do 
#  sample=$(basename $i _aligned_clean_aligned_filtered.bam)
#  echo $sample
#  samtools sort -o ${ALIGNED_WHOLE_DIR}/${sample}_mapped_filtered_sorted.bam $i 
#done

# list the pbmm2 aligned individual files to be merged
# create file if not present
if [ ! -f ${UTILS_DIR}/combined_files.txt ]; then 
  ls /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/pbmm2_align/*Whole*filtered_sorted.bam > ${UTILS_DIR}/combined_files.txt
  ls /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/5_cupcake/5_align/*Targeted*filtered_sorted.bam >> ${UTILS_DIR}/combined_files.txt
fi

# merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y)
cd ${MERGED_CHROM_DIR}
#source activate nanopore
#samtools merge -f WholeTargeted.bam -b ${UTILS_DIR}/combined_files.txt
samtools sort -o WholeTargeted.sorted.bam WholeTargeted.bam
samtools index WholeTargeted.sorted.bam

