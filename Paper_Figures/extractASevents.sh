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


cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/TargetGenes/
find . -type f -name '*final_transcript_classifications.csv' -exec cat {} + > all_final_transcript_classifications.csv
awk '{print $1}'  /lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_2reads2samples_monomultirem.txt > retainedIso.txt
grep -F retainedIso.txt all_final_transcript_classifications.csv > retained_final_transcript_classifications.csv
