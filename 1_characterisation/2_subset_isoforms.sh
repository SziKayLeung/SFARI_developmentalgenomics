#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_subset_isoforms.o
#SBATCH --error=2_subset_isoforms.e

## print start date and time
echo Job started on:
date -u


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation
source $SC_ROOT/SFARI_characterisation.config
source $SC_ROOT/0_source_functions.sh



##-------------------------------------------------------------------------


# output directory 
mkdir -p ${SQ_SUBSET_DIR}

# Conditions to split the files 
cond=(Case Control Case_Control)


##-------------------------------------------------------------------------

source activate nanopore

### Subset by condition and by sample
# subset cases and control by counts
Rscript $SQTABCOUNTS -f ${BRAIN_CLASS} -m ${SAMPLE_META} -c ${BRAIN_ABUNDANCE} -o ${SQ_SUBSET_DIR}
Rscript $SQCOUNT -f ${SQ_SUBSET_DIR}/sqanti3Filtered2_classification.filtered_lite_classification_counts.txt -m ${SAMPLE_META} -o ${SQ_SUBSET_DIR} -d partial

### generate classification.txt, gtf file, and corresponding SQANTI report
sq_dir=$(dirname ${BRAIN_CLASS})
sqname=$(basename ${BRAIN_CLASS})
sqname_prefix=${sqname//"_classification.txt"/}

# for each condition (case, control, case_control)
for c in ${cond[@]}; do 

  echo "Processing $c"
  
  Rscript $SQSUBSET -i ${SQ_SUBSET_DIR}/$c"_ID.txt" -d $sq_dir -s $sqname_prefix -n $c -o ${SQ_SUBSET_DIR} -f 1
  Rscript $SQ_Report ${SQ_SUBSET_DIR}/$c"_classification.txt" ${SQ_SUBSET_DIR}/$c"_junctions.txt"
  
  if [ $c == "Case_Control" ] ; then
    # for Case and control, convert gtf to bed file 
    convert_gtf_bed12 ${SQ_SUBSET_DIR}/$c.gtf 
    source activate nanopore
    
    # colour the bed file by abundance 
    python $ISOCOL --bed ${SQ_SUBSET_DIR}/$c"_sorted.bed12" --a ${SQ_SUBSET_DIR}/$c"_CaseAbundance.csv" --o ${SQ_SUBSET_DIR}/$c"sorted_Case_Coloured.bed12"
    python $ISOCOL --bed ${SQ_SUBSET_DIR}/$c"_sorted.bed12" --a ${SQ_SUBSET_DIR}/$c"_ControlAbundance.csv" --o ${SQ_SUBSET_DIR}/$c"sorted_Control_Coloured.bed12"
  fi

done




