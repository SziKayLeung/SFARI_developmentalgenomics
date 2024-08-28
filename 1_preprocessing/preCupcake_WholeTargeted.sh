#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

# 20/08/2024: run cupcake collapse on merged whole and targeted dataset (timed out)
# 28/08/2024: run cupcake collapse on merged whole and targeted dataset split by chromosome

##-------------------------------------------------------------------------

echo Job started on:
date -u

module load Miniconda2/4.3.21
source activate nanopore

# paths
UTILS_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/0_utilities
MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
MERGED_CUPCAKE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/cupcakeMerged

##-------------------------------------------------------------------------

# merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y)
cd ${MERGED_CHROM_DIR}
samtools merge -f WholeTargeted.bam -b ${UTILS_DIR}/combined_files.txt
samtools index ${MERGED_CHROM_DIR}/WholeTargeted.bam
