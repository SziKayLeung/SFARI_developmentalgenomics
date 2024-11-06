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
#SBATCH --output=1_run_proteogenomics.o
#SBATCH --error=1_run_proteogenomics.e

# 04/12/2023: run on whole dataset
# 04/11/2024: re-run on whole dataset after subsetting from whole+targeted


#-----------------------------------------------------------------------#
## print start date and time
echo Job started on:
date -u

module load Miniconda2
source activate sqanti2_py3
scriptDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/2_characterisation/2_proteinLevel
source ${scriptDir}/0_proteogenomics_functions.sh
source ${scriptDir}/0_proteomics.config

echo "#************************************* Collate and prepare long-read data"
prepare_reference_tables
summarise_longread_data

echo "#************************************* Call open reading frames and classify proteins"
call_orf
refine_calledorf  
classify_protein

echo "#***************All done!****************#"
