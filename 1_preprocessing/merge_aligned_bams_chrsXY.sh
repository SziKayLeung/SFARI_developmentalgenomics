#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=24:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=20G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=v.chundru@exeter.ac.uk # email address


samtools merge -o /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_cleaned_aligned_merged_chrX.bam -b /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_files.txt -R chrX

samtools index /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_cleaned_aligned_merged_chrX.bam

samtools merge -o /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/WholeTargeted_cleaned_aligned_merged_chrX.bam -b /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/combined_files.txt -R chrX

samtools index /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/WholeTargeted_cleaned_aligned_merged_chrX.bam

samtools merge -o /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_cleaned_aligned_merged_chrY.bam -b /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_files.txt -R chrY

samtools index /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_cleaned_aligned_merged_chrY.bam

samtools merge -o /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/WholeTargeted_cleaned_aligned_merged_chrY.bam -b /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/combined_files.txt -R chrY

samtools index /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/WholeTargeted_cleaned_aligned_merged_chrY.bam
