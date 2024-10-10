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
#SBATCH --output=log/demux-%A_%a.o
#SBATCH --error=log/demux-%A_%a.e
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

##-------------------------------------------------------------------------

MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
TCLEAN_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/6_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing


##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
source activate nanopore

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
chrNum=${chromosomes[${SLURM_ARRAY_TASK_ID}]}  

echo $chrNum

if [ -f ${DEMUX_DIR}/WholeTargeted_demux_${chrNum}_fl_count.csv ]; then

  echo WholeTargeted_demux_${chrNum}_fl_count.csv present
  
else

  echo Run demux for $chrNum
  demux_cupcake_collapse.py ${MERGED_CHROM_DIR}/${chrNum}.read_stat.renamed.txt ${DEMUX_DIR}/WholeTargeted_sample_id.csv -o WholeTargeted_demux_${chrNum} --dataset ont -d ${DEMUX_DIR}
  
fi
