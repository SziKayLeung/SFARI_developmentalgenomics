#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=120:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion


##-------------------------------------------------------------------------

MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
TCLEAN_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/6_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing


##-------------------------------------------------------------------------

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${TCLEAN_DIR} -o Targeted -i clean.fa -d ${DEMUX_DIR}
#adapt_cupcake_to_ont.py ${TCLEAN_WHOLE_DIR} -o Whole -i fa -d ${DEMUX_DIR}

# remoe "combined_sorted" string in Targeted_sample_id.csv (due to the input of the fasta files)
cd ${DEMUX_DIR}
awk -F',' '{gsub("_combined_sorted", "", $2); print $1 "," $2}' Targeted_sample_id.csv > Targeted_sample_id_mod.csv

# concatenate whole and targeted sample_id.csv for downstream demux
cat Whole_sample_id.csv Targeted_sample_id_mod.csv  > WholeTargeted_sample_id.csv
