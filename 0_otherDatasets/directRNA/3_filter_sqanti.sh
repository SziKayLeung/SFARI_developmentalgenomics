#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=40:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion

# 29/10/2024: generate fasta, gtf for filtered_sqanti.R

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/6_isoseq_collapse
TCLEAN_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/4_TranscriptClean
SQANTI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/7_sqanti
SQANTI_FINAL_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/9_sqanti_final

module load Miniconda2
source activate nanopore 

#### ------ final filters 

cd ${SQANTI_FINAL_DIR}
outputName=sqantifiltered_monoexonicfiltered_2reads2samples
filteredID=sqantifiltered_monoexonicfiltered_2reads2samples_ID.txt
# fasta
seqtk subseq ${SQANTI_DIR}/SQANTI_merged_all_collapsed_corrected.fasta ${filteredID} > ${outputName}_corrected.fasta
# gtf
grep -wF -f ${filteredID} ${SQANTI_DIR}/SQANTI_merged_all_collapsed_corrected.gtf > ${outputName}.filtered.gtf
# gff
grep -wF -f ${filteredID} ${SQANTI_DIR}/SQANTI_merged_all_collapsed_corrected.gtf.cds.gff > ${outputName}_corrected.gtf.cds.gff
# junctions
grep -wF -f ${filteredID} ${SQANTI_DIR}/SQANTI_merged_all_collapsed_junctions.txt > ${outputName}_junctions.txt

# replace PB ID
python /lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous/replace_PB_id.py \
	-c=${SQANTI_FINAL_DIR}/sqantifiltered_monoexonicfiltered_2reads2samples_classification.txt \
	-g=${SQANTI_FINAL_DIR}/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf \
	-f=${SQANTI_FINAL_DIR}/sqantifiltered_monoexonicfiltered_2reads2samples_corrected.fasta
	
python /lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous/replace_PB_id.py \
	-c=${sqanti}/sqantifiltered_monoexonicfiltered_classification.txt 