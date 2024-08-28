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
#SBATCH --array=0-23%5 # 24 chromosomes
#SBATCH --output=findDuplicatedReads-%A_%a.o
#SBATCH --error=findDuplicatedReads-%A_%a.e

#SBATCH --mem=200G # specify bytes memory to reserve


##-------------------------------------------------------------------------

# directory paths
inputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/WholeTargeted
inputDirStats=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/WholeTargeted/cleaned_merged_collapsed
outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/2_trimmed/targeted/duplicated

# Define the path to the mergedReadStat file from cupcake collapse
clusterReports=(${inputDir}/Targeted*cluster_report.csv)
ReadStats=(${inputDirStats}/WholeTargeted_cleaned_aligned_merged_collapsed_chr*.read_stat.renamed.txt)

##-------------------------------------------------------------------------

wc -l ${clusterReports[@]} >  ${outputDir}/targeted_totalNumber_Reads.csv
sed -i 's#/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/WholeTargeted/##g' ${outputDir}/targeted_totalNumber_Reads.csv
sed -i 's/\.fa_cluster_report\.csv//g' ${outputDir}/targeted_totalNumber_Reads.csv
#for i in ${clusterReports[@]}; do 
#  echo $i
#  sample=$(basename ${i} | cut -d "." -f 1 )
#  sort $i | uniq -cd | awk -F',' '{print $2}' > ${outputDir}/${sample}_duplicated_reads.csv
#done
wc -l ${outputDir}/*_duplicated_reads.csv > ${outputDir}/targeted_totalNumber_duplicatedReads.csv
sed -i 's#/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/2_trimmed/targeted/duplicated/##g' ${outputDir}/targeted_totalNumber_duplicatedReads.csv
sed -i 's/\_duplicated_reads\.csv//g' ${outputDir}/targeted_totalNumber_duplicatedReads.csv

ReadStat=${ReadStats[${SLURM_ARRAY_TASK_ID}]}
chromosome=$(basename ${ReadStat} | cut -d "." -f 1 )
echo $ReadStat
echo $chromosome
cut -f1 ${ReadStat} | sort | uniq -d > ${outputDir}/${chromosome}_duplicated_read_stats.csv 
awk 'NR==FNR{dups[$1]; next} $1 in dups' ${outputDir}/${chromosome}_duplicated_read_stats.csv ${ReadStat} > ${outputDir}/${chromosome}_duplicated_read_stats_final.csv 
