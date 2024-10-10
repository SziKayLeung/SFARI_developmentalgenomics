#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../log/log_Oct2024/3_split_chromosome_for_sqanti.o
#SBATCH --error=../log/log_Oct2024/3_split_chromosome_for_sqanti.e

# 03/09/2024: split chromosomes for downstream sqanti (note chromosome 19 still running)
# 10/10/2024: re-run due to updated transcriptclean paramater to --maxindelLength=10

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
SPLIT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/split

for i in ${ISOSEQ_COLLAPSE_DIR}/*read_stat.txt; do 
  chromosome=$(basename $i .read_stat.txt)
  echo $chromosome
  replaceONTprefix="${chromosome//chr/ONT}"
  #echo $replaceONTprefix
  
  cd ${SPLIT}
  
  # extract the isoform list from each read stat file
  awk -F'\t' '{print $2}' $i | tail -n +2 > WholeTargeted_${chromosome}.id.txt
  
  # replace PB with ONT and the chromosome 
  sed -i "s/PB/${replaceONTprefix}/g" WholeTargeted_${chromosome}.id.txt
  
  # replace gff with the correct gff name
  sed "s/PB/${replaceONTprefix}/g" $ISOSEQ_COLLAPSE_DIR/$chromosome.gff > $ISOSEQ_COLLAPSE_DIR/$chromosome.renamed.gff
  
  # replace read.stat txt with corrected isoform id
  sed "s/PB/${replaceONTprefix}/g" $ISOSEQ_COLLAPSE_DIR/$chromosome.read_stat.txt > $ISOSEQ_COLLAPSE_DIR/$chromosome.read_stat.renamed.txt
  
  # split chunk with 1000000 lines
  split WholeTargeted_${chromosome}.id.txt WholeTargeted_${chromosome}_chunk -l1000000 --additional-suffix=.txt -d

done

ls *chunk* > WholeTargeted_chunks.txt