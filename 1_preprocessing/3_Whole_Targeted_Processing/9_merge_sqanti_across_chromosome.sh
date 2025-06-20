#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --mem=200G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

MERGED_SQANTI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax_merged
SQANTI_FINAL_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous

module load Miniconda2
source activate nanopore 

# merge classification files 
echo "******************** classification files"
ls ${MERGED_SQANTI_DIR}/*RulesFilter_result_classification.txt*
classFiles=($(ls ${MERGED_SQANTI_DIR}/*RulesFilter_result_classification.txt*))
# Concatenate the header from the first file and content from all files, skipping headers in subsequent files
{ head -n 1 "${classFiles[0]}"; tail -n +2 -q "${classFiles[@]}"; } > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr_RulesFilter_result_classification.txt

# merge filtered.gtf across multiple chromosomes
echo "******************** filtered gtf"
ls ${MERGED_SQANTI_DIR}/*filtered.gtf
filtered_gtf=($(ls ${MERGED_SQANTI_DIR}/*filtered.gtf))
cat ${filtered_gtf[@]} > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr.filtered.gtf

# merge cds gff across multiple chromosomes
echo "******************** filtered gff"
ls ${MERGED_SQANTI_DIR}/*cds.gff
cds_gff=($(ls ${MERGED_SQANTI_DIR}/*cds.gff))
cat ${cds_gff[@]} > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr_corrected.gtf.cds.gff

# merge corrected.fasta across multiple chromosomes
echo "******************** fasta"
ls ${MERGED_SQANTI_DIR}/*corrected.fasta
corrected_fasta=($(ls ${MERGED_SQANTI_DIR}/*corrected.fasta))
cat ${corrected_fasta[@]} > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr_corrected.fasta

# merge junctions.txt 
echo "******************** junctions"
ls ${MERGED_SQANTI_DIR}/*junctions*
junctions=($(ls ${MERGED_SQANTI_DIR}/*junctions*))
# Concatenate the header from the first file and content from all files, skipping headers in subsequent files
{ head -n 1 "${junctions[0]}"; tail -n +2 -q "${junctions[@]}"; } > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr_junctions.txt

# merge filtering reasons 
echo "******************** filtering reasons"
ls ${MERGED_SQANTI_DIR}/*filtering_reasons.txt*
filteringReasons=($(ls ${MERGED_SQANTI_DIR}/*filtering_reasons.txt*))
# Concatenate the header from the first file and content from all files, skipping headers in subsequent files
{ head -n 1 "${filteringReasons[0]}"; tail -n +2 -q "${filteringReasons[@]}"; } > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsedAllChr_filtering_reasons.txt


#### ------ final filters 

cd ${SQANTI_FINAL_DIR}
outputName=sqantifiltered_monoexonicfiltered_2reads2samples
filteredID=sqantifiltered_monoexonic_2reads2samplesfiltered_ID.txt
# fasta
seqtk subseq WholeTargeted_collapsedAllChr_corrected.fasta ${filteredID} > ${outputName}_corrected.fasta
# gtf
grep -wF -f ${filteredID} WholeTargeted_collapsedAllChr.filtered.gtf > ${outputName}.filtered.gtf
# gff
grep -wF -f ${filteredID} WholeTargeted_collapsedAllChr_corrected.gtf.cds.gff > ${outputName}_corrected.gtf.cds.gff
# junctions
grep -wF -f ${filteredID} WholeTargeted_collapsedAllChr_junctions.txt > ${outputName}_junctions.txt

