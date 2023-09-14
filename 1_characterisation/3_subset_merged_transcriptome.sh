#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../log/3_subset_merged_transcriptome.o
#SBATCH --error=../log/3_subset_merged_transcriptome.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 

wholeTargetedDir=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/
wholeTargetedGff=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/merged_sorted.gff
WK_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/

source activate sqanti2_py3
#Rscript ${SC_ROOT}/3_subset_merged_transcriptome.R
subset_fasta_gtf.py ${wholeTargetedGff} --gtf -i ${wholeTargetedDir}/filteredIds.csv -o merged_sorted_filtered2reads2samples -d ${WK_DIR}