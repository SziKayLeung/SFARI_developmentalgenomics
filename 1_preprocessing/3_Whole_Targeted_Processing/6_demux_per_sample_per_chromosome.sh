#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --array 0-23 # 24 chromsomes, 22 autosomal, X and Y
#SBATCH --output=../log/log_Oct2024/6_demux_persample_perchromosome-%A_%a.o
#SBATCH --error=../log/log_Oct2024/6_demux_persample_perchromosome-%A_%a.e
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
source activate nanopore

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
chrNum=${chromosomes[${SLURM_ARRAY_TASK_ID}]}  

for sample_id in  ${DEMUX_DIR}/1_ID/*unique_sample_id.csv*; do 
	filename=$(basename $sample_id)
	echo $filename
	prefix="${filename%%_*}"
	echo $prefix
	demux_cupcake_collapse.py ${DEMUX_DIR}/2_finalised_readstat/${chrNum}.read_stat.renamed.min2FL.sqantiMonoexonicFiltered.txt ${sample_id} -o ${prefix}_demux_${chrNum} --dataset ont -d ${DEMUX_DIR}/3_demux
done 

mkdir ${DEMUX_DIR}/4_demux_merged
demux_post_samplechromosome.py -i=${DEMUX_DIR}/3_demux --chr=${chrNum} -o_name=WholeTargeted -o_dir=${DEMUX_DIR}/4_demux_merged
