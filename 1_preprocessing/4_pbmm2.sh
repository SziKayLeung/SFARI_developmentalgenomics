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
#SBATCH --array=0-39%5
#SBATCH --output=4_pbmm2_filter-%A_%a.o
#SBATCH --error=4_pbmm2_filter-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
source activate nanopore

# set up batch
WKD_ROOT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome
#
RAW_FASTQ=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
export TCLEAN=/lustre/projects/Research_Project-MRC148213/lsl693/software/TranscriptClean/TranscriptClean.py

porechop_fastq_files=($(ls ${RAW_FASTQ}/*fastq.gz))
SamplePath=${porechop_fastq_files[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${SamplePath} .fastq.gz)
echo "Processing ${sample}"

source activate isoseq3

# pbmm2 align
cd ${WKD_ROOT}/5_cupcake/5_align
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${WKD_ROOT}/4_tclean/${sample}_combined_sorted.sam/${sample}_combined_sorted.sam_clean.fa ${sample}_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log
samtools view -F 4032 ${sample}_mapped.bam -b -o ${sample}_mapped_filtered.bam
samtools sort -o ${sample}_mapped_filtered_sorted.bam ${sample}_mapped_filtered.bam 
samtools index ${sample}_mapped_filtered_sorted.bam