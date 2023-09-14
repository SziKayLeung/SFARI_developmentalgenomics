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


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation
source $SC_ROOT/SFARI_characterisation.config
source $SC_ROOT/01_source_functions.sh


cp ${BRAIN_SQ_WHOLE_GTF} ${CHAR_DIR}
convert_gtf_bed12 ${CHAR_DIR}/sqanti3Filtered2_classification.filtered_lite.gtf
cd ${CHAR_DIR}
grep ENSG00000186868 sqanti3Filtered2_classification.filtered_lite.gtf > MAPT_sqanti3Filtered_classification.filtered_lite.gtf

