#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=60:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --array 0-23 # 24 chromsomes, 22 autosomal, X and Y
#SBATCH --output=../log/log_Oct2024/6_demux-%A_%a.o
#SBATCH --error=../log/log_Oct2024/6_demux-%A_%a.e
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
TCLEAN_TARGETED_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/B_Targeted/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing


##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
source activate nanopore

# set batch array (run per chromosome)
renamedFiles=(${ISOSEQ_COLLAPSE_DIR}/split/*renamed*)
renamedFile=${renamedFiles[$SLURM_ARRAY_TASK_ID]}
fileName=$(basename $renamedFile)

echo $renamedFile
demux_cupcake_collapse.py $renamedFile ${DEMUX_DIR}/WholeTargeted_sample_id.csv -o WholeTargeted_demux_fileName --dataset ont -d ${DEMUX_DIR}
