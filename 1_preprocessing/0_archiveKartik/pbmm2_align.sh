#!/bin/bash

#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=120:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16# specify number of processors per node
#SBATCH --mem=250G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=v.chundru@exeter.ac.uk # email address
#SBATCH --array 44 
#27,41 #1-43,45-87


i=`sed -n "${SLURM_ARRAY_TASK_ID}p;d" /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/bamfiles_to_clean.txt`

j=`echo ${i} | tr '/' ' ' | awk '{print $8}'`

pbmm2 align --preset ISOSEQ --sort /lustre/home/vc362/resources/Homo_sapiens_assembly38.fasta ${i}_clean.fa /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/${j}_clean_aligned.bam --log-level TRACE --log-file /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/${j}.log --num-threads 12 --sort-threads 4

samtools view -F 4032 /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/${j}_clean_aligned.bam -b -o /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford/pbmm2_align/${j}_clean_aligned_filtered.bam
