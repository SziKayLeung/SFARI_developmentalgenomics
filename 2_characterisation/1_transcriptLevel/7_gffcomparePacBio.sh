#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

# 13/02/2024: overlap between PacBio whole cortex dataset 
module load Miniconda2/4.3.21

sqantiDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti
sqanti_gtf=WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_Whole_2reads2samples_nomonointergenic.gtf

outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/14_OverlapPacBio
PacBioDir=/lustre/projects/Research_Project-MRC148213/lsl693/PacBioPaper/SQANTI2/HumanCTX
PacBio_gtf=HumanCTX.collapsed_classification.filtered_lite.gtf

## ----- run gffcompare ---- 

source activate sqanti2_py3 
cd ${outputDir}
cp ${sqantiDir}/${sqanti_gtf} .
cp ${PacBioDir}/${PacBio_gtf} .
PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r ${sqanti_gtf} ${PacBio_gtf} -o sfariPacBio
gffcompare -r ${PacBio_gtf} ${sqanti_gtf} -o PacBiosfari
