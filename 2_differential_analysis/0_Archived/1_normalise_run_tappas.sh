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
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/2_differential_analysis
source $SC_ROOT/sfari_differential.config


## ---------------- A) Fetal vs Adult Brain

source activate nanopore 
# extract expression matrix from talon abundance file
Rscript $GENERALFUNC/5_TappAS_Differential/talon2tappas_expression.R -e ${EXP} -o $TAPPAS_INPUT_DIR}/A_FetalvsAdult -n ${NAME}"_talon_expression"

# normalise
Rscript $GENERALFUNC/5_TappAS_Differential/normalise.R -e $TAPPAS_INPUT_DIR}/A_FetalvsAdult/${NAME}"_talon_expression.txt" -f ${TAPPAS_FACTOR} -s ${SQ_CLASS} -o $TAPPAS_INPUT_DIR}/A_FetalvsAdult -n ${NAME} 

# prepare tappAS files 
cp $SQ_GFF3 $TAPPAS_INPUT_DIR}/A_FetalvsAdult

# run on knight 
# perform differential expression analysis
scp sl693@login.isca.ex.ac.uk:${TAPPAS_INPUT_DIR}/A_FetalvsAdult/* /mnt/data1/Szi/SFARI/A_FetalvsAdult

# transfer tappAS results back to ISCA 
# normalised externally
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.1523097622.tappas/* ${TAPPAS_OUTPUT_DIR}/A_FetalvsAdult

# normalised internally using tappAS java 
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.0842801692.tappas/* ${TAPPAS_OUTPUT_DIR}/A_FetalvsAdult


## ---------------- B) Fetal vs Adult Brain (Sex differences)

# copy sqanti annotation file 
cp $SQ_GFF3 $TAPPAS_INPUT_DIR}/B_FetalvsAdultvsSex

# data wrangle and copy the expression file
Rscript $GENERALFUNC/5_TappAS_Differential/talon2tappas_expression.R -e ${EXP} -o ${TAPPAS_INPUT_DIR}/B_FetalvsAdultvsSex -n ${NAME}"_talon_expression"

# run on knight 
# perform differential expression analysis
scp -r sl693@login.isca.ex.ac.uk:${TAPPAS_INPUT_DIR}/B_FetalvsAdultvsSex/* /mnt/data1/Szi/SFARI/B_FetalvsAdultvsSex

# transfer tappAS results back to ISCA 
# normalised externally
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.01592099705.tappas/* ${TAPPAS_OUTPUT_DIR}/B_FetalAdultvsSex
