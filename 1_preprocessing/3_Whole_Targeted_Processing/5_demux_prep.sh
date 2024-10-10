#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=40:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=../log/log_Oct2024/5_demux_prep.o
#SBATCH --error=../log/log_Oct2024/5_demux_prep.e


##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
TCLEAN_TARGETED_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/B_Targeted/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing


##-------------------------------------------------------------------------

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${TCLEAN_TARGETED_DIR} -o Targeted -i clean.fa -d ${DEMUX_DIR}
adapt_cupcake_to_ont.py ${TCLEAN_WHOLE_DIR} -o Whole -i fa -d ${DEMUX_DIR}

# remoe "combined_sorted" string in Targeted_sample_id.csv (due to the input of the fasta files)
cd ${DEMUX_DIR}
awk -F',' '{gsub("_combined_sorted", "", $2); print $1 "," $2}' Targeted_sample_id.csv > Targeted_sample_id_mod.csv

# concatenate whole and targeted sample_id.csv for downstream demux
cat Whole_sample_id.csv Targeted_sample_id_mod.csv  > WholeTargeted_sample_id.csv
