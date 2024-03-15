#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=1b_extract_monoexonic_transcripts.e
#SBATCH --output=1b_extract_monoexonic_transcripts.o

# 01/02/2024: Extract monoexonic transcripts from reference human genocode annotations


##-------------------------------------------------------------------------

LOGEN=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
REF=/lustre/projects/Research_Project-MRC148213/lsl693/reference/human
OUTDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/0_utils
export PATH=$PATH:${LOGEN}/miscellaneous/


##-------------------------------------------------------------------------

module load Miniconda2
source activate nanopore
extract_exon_num.py ${REF}/gencode.v40.annotation.gtf -d ${OUTDIR}