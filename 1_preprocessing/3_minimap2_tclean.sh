#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-39%5
#SBATCH --output=3_minimap2_tclean-%A_%a.o
#SBATCH --error=3_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
source activate nanopore

# set up batch
WKD_ROOT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome
RAW_FASTQ=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
export TCLEAN=/lustre/projects/Research_Project-MRC148213/lsl693/software/TranscriptClean/TranscriptClean.py

porechop_fastq_files=($(ls ${RAW_FASTQ}/*fastq.gz))
SamplePath=${porechop_fastq_files[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${SamplePath} .fastq.gz)
echo "Processing ${sample}"

# 6) run_minimap2 <input_fasta> <output_dir>
# Aim: Align reads from trimming, filtering to genome of interest using Minimap2
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_minimap2(){
  
  source activate nanopore
  
  name=$(basename $1 .fasta)
  echo "Aligning ${name} using Minimap2"
  
  minimap2 -t 46 -ax splice ${GENOME_FASTA} $1 > $2/${name}.sam 2> $2/${name}_minimap2.log
  samtools sort -O SAM $2/${name}.sam > $2/${name}_sorted.sam
  
  source deactivate
}


# run_transcriptclean <input_sam> <output_dir>
run_transcriptclean(){
  
  source activate sqanti2_py3
  
  name=$(basename $1 _merged_combined_sorted.sam)
  echo "TranscriptClean ${name}"
  
  cd $2; mkdir -p ${name}
  cd $2/${name}
  python ${TCLEAN} --sam $1 --genome ${GENOME_FASTA} --outprefix $2/${name}/${name} --tmpDir $2/${name}/${name}_tmp
}


##-------------------------------------------------------------------------

# map combined fasta to reference genome
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/${sample}_combined.fasta ${WKD_ROOT}/3_minimap

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/3_minimap/${sample}_combined_sorted.sam ${WKD_ROOT}/4_tclean

source activate isoseq3

# pbmm2 align
cd ${WKD_ROOT}/5_cupcake/5_align
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${WKD_ROOT}/4_tclean/${sample}_combined_sorted.sam/${sample}_combined_sorted.sam_clean.fa ${sample}_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log