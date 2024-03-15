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
#SBATCH --error=HervGffCompare.e
#SBATCH --output=HervGffCompare.o

# 13/01/2024: gffcompare on unfiltered collapsed dataset and HERV annotations

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
sqantiDir=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3
utilsDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/0_utils
unfilteredGtf=${sqantiDir}/unfiltered/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected.gtf 
hervGtf=${utilsDir}/HERV.gtf
outputDir=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/8_characterisation/0_geneSpecific/HERV

## ----- run gffcompare ---- 

#source activate sqanti2_py3 
cd ${outputDir}
#PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
#gffcompare -r $unfilteredGtf $hervGtf -o HERVUnfiltered

#cd $utilsDir
#source activate sqanti2_py3
#gtfToGenePred HERV.gtf HERV.genePred
#genePredToBed HERV.genePred > HERV.bed12

#cd ${outputDir}
#bedtools intersect -wa -wb -a ${unfilteredGtf} -b ${utilsDir}/HERV.bed12 > overlap.bed

output_file="overlap_isoform.bed"
while IFS= read -r line; do
  # Extract ONTXX value using grep and sed
  ont_value=$(echo "$line" | grep -o -P 'ONT\d+_\d+_\d+' | sed 's/"//g')

  # Print the extracted ONTXX value to the output file
  echo "$ont_value" >> "$output_file"
done < overlap.bed

grep -f overlap_isoform.bed ${unfilteredGtf}  > overlap.gtf
awk -F'\t' 'BEGIN {OFS=FS} {gsub(/"/, "", $9); print}' overlap.gtf > overlapFinal.gtf