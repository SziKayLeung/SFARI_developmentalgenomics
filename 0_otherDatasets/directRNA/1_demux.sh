#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=40:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion

# 29/10/2024: demux output for direct RNA sequencing dataset

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/6_isoseq_collapse
TCLEAN_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/4_TranscriptClean
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/8_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing


##-------------------------------------------------------------------------

module load Miniconda2
source activate nanopore

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${TCLEAN_DIR} -o dRNA -i clean.fa -d ${DEMUX_DIR}

# demux script
demux_cupcake_collapse.py ${ISOSEQ_COLLAPSE_DIR}/merged_all_collapsed.read_stat.txt ${DEMUX_DIR}/dRNA_sample_id.csv -o merged --dataset ont -d ${DEMUX_DIR}