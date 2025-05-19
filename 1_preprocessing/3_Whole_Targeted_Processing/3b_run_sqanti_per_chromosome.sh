#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../log/log_Oct2024/3b_run_sqanti_per_chromosome-%A_%a.o
#SBATCH --error=../log/log_Oct2024/3b_run_sqanti_per_chromosome-%A_%a.e

# 30/10/2024: default sqanti json file

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
SQANTI_RELAXED=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax_merged
SQANTI_DEFAULT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_default_merged
SQANTI_FINAL_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final

## Software 
export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export CUPCAKE=${SOFTDIR}/cDNA_Cupcake
export SEQUENCE=$CUPCAKE/sequence
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export SQANTI3_DATA_DIR=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data
export CAGE_PEAK=${SQANTI3_DATA_DIR}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=${SQANTI3_DATA_DIR}/polyA_motifs/mouse_and_human.polyA_motif.txt
export refAnno=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/gencode.v38.annotation.gtf
export refFile=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/Hg38.fa
export SQANTI_JSON=$SQANTI3_DIR/utilities/filter/filter_default.json

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
 
echo "Run Sqanti on $chrNum"
source activate sqanti2_py3
  
## run sqanti filter
for chrNum in ${chromosomes[@]}; do 
  
  echo "Run Sqanti filter $chrNum"
  cd ${SQANTI_DEFAULT}
  cp ${SQANTI_RELAXED}/WholeTargeted_collapsed${chrNum}_corrected.gtf . 
  cp ${SQANTI_RELAXED}/WholeTargeted_collapsed${chrNum}_classification.txt .
  
  python ${SQANTI3_DIR}/sqanti3_filter.py rules \
  	  --output WholeTargeted_collapsed${chrNum} \
  	  --skip_report --gtf WholeTargeted_collapsed${chrNum}_corrected.gtf \
  	  --json_filter ${SQANTI_JSON} WholeTargeted_collapsed${chrNum}_classification.txt
     
done


echo "******************** filtered gtf"
ls ${SQANTI_DEFAULT}/*filtered.gtf
filtered_gtf=($(ls ${SQANTI_DEFAULT}/*filtered.gtf))
cat ${filtered_gtf[@]} > ${SQANTI_FINAL_DIR}/WholeTargeted_collapsed_AllChrDefaultSqanti.filtered.gtf
