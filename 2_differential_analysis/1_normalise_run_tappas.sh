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


## ---------------------------

source activate nanopore 
# extract expression matrix from talon abundance file
Rscript $GENERALFUNC/5_TappAS_Differential/talon2tappas_expression.R -e ${EXP} -o ${TAPPAS_DIR} -n ${NAME}"_talon_expression"

# normalise
Rscript $GENERALFUNC/5_TappAS_Differential/normalise.R -e ${TAPPAS_DIR}/${NAME}"_talon_expression.txt" -f ${TAPPAS_FACTOR} -s ${SQ_CLASS} -o ${TAPPAS_DIR} -n ${NAME} 

# prepare tappAS files 
cp $SQ_GFF3 ${TAPPAS_DIR}

# run on knight 
scp sl693@login.isca.ex.ac.uk:/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/9_tappAS_SK/* /mnt/data1/Szi/SFARI