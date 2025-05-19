#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
filteredID=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_output/sqantifiltered_monoexonicfiltered_ID.txt
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux

for i in ${ISOSEQ_COLLAPSE_DIR}/*read_stat.txt; do 
  	chromosome=$(basename $i .read_stat.txt)
  	echo ${chromosome}
	  # Get the header from the source file
	  header=$(head -n 1 ${ISOSEQ_COLLAPSE_DIR}/${chromosome}.read_stat.renamed.min2FL.txt)

  	# Print the header to the output file
	  echo "${header}" > ${DEMUX_DIR}/2_finalised_readstat/${chromosome}.read_stat.renamed.min2FL.sqantiMonoexonicFiltered.txt
  
	  # Append the filtered content below the header
	  grep -Ff ${filteredID} ${ISOSEQ_COLLAPSE_DIR}/${chromosome}.read_stat.renamed.min2FL.txt >> ${DEMUX_DIR}/2_finalised_readstat/${chromosome}.read_stat.renamed.min2FL.sqantiMonoexonicFiltered.txt
done 