#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# 7/11/2024: merge output from mass-spec alignment

#-----------------------------------------------------------------------

module load Miniconda2
source activate sqanti2_py3
source /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/2_characterisation/2_proteinLevel/0_mass_spec.config

cd ${Part2WKD}
sfaripeptides=$(find -L . -name "Sfari_peptides.bed12")
cat $sfaripeptides > merged_sfari_peptides.bed12

find . -type f -name *_novel_peptides.gtf* -exec cat {} + > AllSfari_novelpeptides.gtf
find . -type f -name '*novel_peptides_seq.bed12*' -exec cat {} + > Allsfari_novel_peptides.bed12
sort Allsfari_novel_peptides.bed12 | uniq > Allsfari_novel_peptides_unique.bed12
sort AllSfari_novelpeptides.gtf | uniq > AllSfari_novelpeptides_unique.gtf