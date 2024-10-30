#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=40:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=../log/log_Oct2024/5_demux_prep.o
#SBATCH --error=../log/log_Oct2024/5_demux_prep.e


##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
TCLEAN_TARGETED_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/B_Targeted/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux
ALIGNED_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/5_isoseq/pbmm2_align
ALIGNED_TARGETED_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/B_Targeted/5_cupcake/5_align_filter

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

mkdir -p ${DEMUX_DIR}/1_ID

##-------------------------------------------------------------------------

module load Miniconda2 

source activate nanopore 

for i in ${ALIGNED_WHOLE_DIR}/*filtered_sorted.bam; do 
  echo $i
  filename=$(basename $i _mapped_filtered_sorted.bam)
  samtools bam2fq $i | seqtk seq -A - > ${ALIGNED_WHOLE_DIR}/${filename}_mapped_filtered_sorted.fa
done 

for i in ${ALIGNED_TARGETED_DIR}/*filtered_sorted.bam; do 
  echo $i
  filename=$(basename $i _mapped_filtered_sorted.bam)
  samtools bam2fq $i | seqtk seq -A - > ${ALIGNED_TARGETED_DIR}/${filename}_mapped_filtered_sorted.fa
done 


# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${ALIGNED_TARGETED_DIR} -o Targeted -i mapped_filtered_sorted.fa -d ${DEMUX_DIR}/1_ID
adapt_cupcake_to_ont.py ${ALIGNED_WHOLE_DIR} -o Whole -i mapped_filtered_sorted.fa -d ${DEMUX_DIR}/1_ID
