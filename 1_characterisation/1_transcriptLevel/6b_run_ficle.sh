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


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2
source activate ficle
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/FICLE
export PATH=$PATH:${FICLE_ROOT}


##-------------------------------------------------------------------------

# directory paths
inputDir=/lustre/home/vc362/lustre_project/ficle/input_files/
outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/
classFile=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_2reads2samples_monomultirem.txt
  
# run ficle
gene=HSPA1B
ficle.py --genename=${gene} \
  --reference=${inputDir}/${gene}_gencode.gtf \
  --input_gtf=${inputDir}/${gene}.gtf  \
  --input_class=${classFile}  \
  --output_dir=${outputDir}
  

input_file="/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/remadeOutput/TargetGenes/missingGenes.txt"
# Loop through each row in the file
while IFS=, read -r col1; do
  # Process each column or row as needed
  gene=$col1
  echo $gene
  ficle.py --genename=${gene}  --reference=${inputDir}/${gene}_gencode.gtf  --input_gtf=${inputDir}/${gene}.gtf    --input_class=${classFile}    --output_dir=${outputDir}
done < "$input_file"

cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/15_ficle/remadeOutput/TargetGenes/
find . -type f -name '*final_transcript_classifications.csv' -exec cat {} + > all_final_transcript_classifications.csv
awk '{print $1}'  /lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_2reads2samples_monomultirem.txt > retainedIso.txt
grep -F retainedIso.txt all_final_transcript_classifications.csv > retained_final_transcript_classifications.csv
