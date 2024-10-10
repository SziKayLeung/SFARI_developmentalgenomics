#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array 0-23 # 24 chromsomes, 22 autosomal, X and Y
#SBATCH --output=cupcake_WholeTargeted-%A_%a.o
#SBATCH --error=cupcake_WholeTargeted-%A_%a.e

# 20/08/2024: run cupcake collapse on merged whole and targeted dataset (timed out)
# 28/08/2024: run cupcake collapse on merged whole and targeted dataset split by chromosome
# 10/08/2024: re-run cupcake collapse due to wrong sample


##-------------------------------------------------------------------------

echo Job started on:
date -u

module load Miniconda2/4.3.21

# paths
UTILS_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/0_utilities
MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
MERGED_CUPCAKE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/cupcakeMerged

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
chrNum=${chromosomes[${SLURM_ARRAY_TASK_ID}]}  


##-------------------------------------------------------------------------

# merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y)
cd ${MERGED_CHROM_DIR}

# isoseq collapse the merged file by chromosome
samtools view -b ${MERGED_CHROM_DIR}/WholeTargeted.sorted.bam $chrNum > ${MERGED_CHROM_DIR}/${chrNum}.bam

source activate isoseq3
isoseq3 collapse ${MERGED_CHROM_DIR}/${chrNum}.bam ${chrNum}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16
#isoseq3 collapse ${MERGED_CHROM_DIR}/WholeTargeted.bam WholeTargeted.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16