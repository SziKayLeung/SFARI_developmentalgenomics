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
#SBATCH --array=0-46%5 # 47 samples
#SBATCH --output=rarefaction_whole-%A_%a.o
#SBATCH --error=rarefaction_whole-%A_%a.e
#SBATCH --job-name=whole_transcriptome_rarefaction

# 20/08/2024: rarefaction curves on whole dataset (47 samples)

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
RAREFACTION=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_rarefactionSKL
CUPCAKE=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/cupcakeIndividual
CUPCAKE_WHOLE=($(ls $CUPCAKE/*Whole*.read_stat.txt*))
SamplePath=${CUPCAKE_WHOLE[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${SamplePath} .read_stat.txt)

module load Miniconda2/4.3.21
source activate sqanti2_py3

# make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <working_directory> <sqanti_class_file>
make_file_for_rarefaction(){
  
  CUPCAKE_ANNOTATION=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/annotation
  
  cd $3
  echo "Working with $1"
  prefix=$1
  # make_file_for_subsampling_from_collapsed.py <sample_name_prefix>.input.file <sample_name_prefix>.output.file <sample_name_prefix>.classification.txt
  python $CUPCAKE_ANNOTATION/make_file_for_subsampling_from_collapsed.py -i $2/$prefix -o $1.subsampling &>> $1.makefile.log
  python $CUPCAKE_ANNOTATION/subsample.py --by pbgene --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_pbgene.min_fl_2.txt
  python $CUPCAKE_ANNOTATION/subsample.py --by pbid --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_pbid.min_fl_2.txt
  
  # if running with sqanti_class file
  #python $CUPCAKE_ANNOTATION/subsample.py --by refgene --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refgene.min_fl_2.txt 
  #python $CUPCAKE_ANNOTATION/subsample.py --by refisoform --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refisoform.min_fl_2.txt
  
}

make_file_for_rarefaction $sample $CUPCAKE $RAREFACTION   