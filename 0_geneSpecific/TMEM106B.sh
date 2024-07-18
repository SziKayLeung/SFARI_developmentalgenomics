#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address 
#SBATCH --output=TMEM106B.o
#SBATCH --error=TMEM106B.e

# 24.05.2024: Extracting TMEM106B raw reads from Whole&Targeted dataset

##-------------------------------------------------------------------------

module load Miniconda2
source activate nanopore

## directory paths
sfari=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI
trimmedDIR=${sfari}/2_trimmed
TmemDIR=${sfari}/0_colloborations/TMEM106B

## import files
# TMEM106B_isoforms.txt generated from TMEM106B.R
IsoformList=${TmemDIR}/TMEM106B_isoforms.txt
# read stat file generated from Iso-Seq cupcake collapse
ReadStat=${sfari}/5_isoseq/WholeTargeted/WholeTargeted_cleaned_aligned_merged_collapsed_chr7.read_stat.renamed.txt


##-------------------------------------------------------------------------

cd ${TmemDIR}
echo "Grep reads"
grep -w -f TMEM106B_isoforms.txt ${ReadStat} > TMEM106B_WholeTargeted.read_stat.txt
cat TMEM106B_WholeTargeted.read_stat.txt | awk '{print $1}' > TMEM106B_WholeTargeted.rawReadID.txt

echo "Concatenate fastq files"
ls ${sfari}/2_trimmed/*q*
cat ${sfari}/2_trimmed/*q* > ./Merged_trimmed.fastq

echo "Extract raw reads"
seqtk subseq Merged_trimmed.fastq TMEM106B_WholeTargeted.rawReadID.txt > TMEM106B_WholeTargeted_rawRead.fastq

