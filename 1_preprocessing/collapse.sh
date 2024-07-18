#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=120:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=v.chundru@exeter.ac.uk # email address
#SBATCH --array 1-22

#isoseq3 collapse /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Targeted_cleaned_aligned_merged_chr${SLURM_ARRAY_TASK_ID}.bam /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Targeted_cleaned_aligned_merged_collapsed_chr${SLURM_ARRAY_TASK_ID}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95

#isoseq3 collapse /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/Whole_cleaned_aligned_merged_chr${SLURM_ARRAY_TASK_ID}.bam /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_chr${SLURM_ARRAY_TASK_ID}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16

isoseq3 collapse /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/merged/WholeTargeted_cleaned_aligned_merged_chr${SLURM_ARRAY_TASK_ID}.bam /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_chr${SLURM_ARRAY_TASK_ID}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16
