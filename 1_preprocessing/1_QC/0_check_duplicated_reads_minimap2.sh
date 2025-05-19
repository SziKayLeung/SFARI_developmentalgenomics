#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion


##-------------------------------------------------------------------------

module load Miniconda2
source activate nanopore

ALIGNED_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/B_Targeted/3_minimap
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

##-------------------------------------------------------------------------

cd ${ALIGNED_DIR}
mkdir -p sequences

# fasta sequences
for i in *sorted.sam*; do 
	echo ${i}
	sample=$(basename ${i} _combined_sorted.sam)
	echo ${sample}
	samtools view -S -b $i | samtools bam2fq - | seqtk seq -A > ${sample}.fa
done 

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${ALIGNED_DIR} -o Targeted_aligned -i fa -d ${ALIGNED_DIR}/sequences

# find duplicated reads
awk '{print $1}' Targeted_aligned_id.csv | sort | uniq -d > duplicated_reads.csv
