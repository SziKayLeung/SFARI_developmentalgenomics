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

# normalise brain vs pancreas expression file (too large to run it on tappAS)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis
source $SC_ROOT/sfari_differential.config


source activate nanopore 
Rscript $GENERALFUNC/5_TappAS_Differential/normalise.R \
-e $TAPPAS_INPUT_DIR/C_BrainvsPancreas/merged_brain_pancreas_fetal_targeted_expression.txt \
-f ${BRAIN_PANCREAS_TAPPAS_FACTOR} \
-s ${BRAIN_PANCREAS_CLASS} \
-i ${BRAIN_PANCREAS_TARGETED_ISO} \
-o $TAPPAS_INPUT_DIR/C_BrainvsPancreas \
-n merged_brain_pancreas_fetal_targeted