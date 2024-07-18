#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_subset_isoforms.o
#SBATCH --error=2_subset_isoforms.e

## Subset control and case samples (n = 36) from Brain ONT whole transcriptome data (prenatal, postnatal)
## Subset control and case samples (n = 35) from Brain ONT targeted data (prenatal, postnatal)

## print start date and time
echo Job started on:
date -u


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/1_transcriptLevel
source $SC_ROOT/SFARI_characterisation.config
source $SC_ROOT/01_source_functions.sh


##-------------------------------------------------------------------------


# output directory 
mkdir -p ${BRAIN_SQ_SUBSET_DIR} ${BRAIN_SQ_SUBSET_WHOLE_DIR} ${BRAIN_SQ_SUBSET_TARGETED_DIR}


##-------------------------------------------------------------------------

# subset_case_control <class_file> <meta_file> <abundance_file> <output_dir> <cpat_file> <cpat_noORF> <control_variable> <case_variable>
subset_case_control ${BRAIN_SQ_WHOLE_CLASS} ${BRAIN_WHOLE_META} ${BRAIN_TALON_WHOLE_ABUNDANCE} ${BRAIN_SQ_SUBSET_WHOLE_DIR} ${BRAIN_WHOLE_CPAT} ${BRAIN_WHOLE_CPAT_NOORF} Prenatal Postnatal 
subset_case_control ${BRAIN_SQ_TARGETED_CLASS} ${BRAIN_TARGETED_META} ${BRAIN_TALON_TARGETED_ABUNDANCE} ${BRAIN_SQ_SUBSET_TARGETED_DIR} ${BRAIN_TARGETED_CPAT} ${BRAIN_TARGETED_CPAT_NOORF} Prenatal Postnatal 