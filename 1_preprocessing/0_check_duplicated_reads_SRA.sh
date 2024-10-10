#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=0:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion


##-------------------------------------------------------------------------

SRA=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined
SRA_ID=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combinedID
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

##-------------------------------------------------------------------------

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${SRA} -o Targeted_sample -i fasta -d ${SRA_ID}

# find duplicated reads
awk '{print $1}' Targeted_sample_id.csv | sort | uniq -d > duplicated_reads.csv