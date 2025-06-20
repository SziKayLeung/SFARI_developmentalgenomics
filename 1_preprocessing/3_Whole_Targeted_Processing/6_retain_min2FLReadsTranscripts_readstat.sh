#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

module load Miniconda2
source activate nanopore

for i in ${ISOSEQ_COLLAPSE_DIR}/*read_stat.txt; do 
  chromosome=$(basename $i .read_stat.txt)
  echo ${chromosome}
  python /lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing/filter_cupcake_collapse_stat.py --fl_read=2 ${ISOSEQ_COLLAPSE_DIR}/${chromosome}.read_stat.renamed.txt ${ISOSEQ_COLLAPSE_DIR}/${chromosome}.abundance.renamed.txt
done 
 
