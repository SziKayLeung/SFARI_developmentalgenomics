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

cd ${SQANTI_FINAL_DIR}

outputName=sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron
filteredID=sqantifiltered_monoexonic_2reads2samplesfiltered_whole_intergenicGenicIntron_ID.txt
grep -wF -f ${filteredID} sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf > ${outputName}.filtered.gtf

source activate nanopore 
seqtk subseq WholeTargeted_collapsedAllChr_corrected.fasta ${filteredID} > ${outputName}.fasta