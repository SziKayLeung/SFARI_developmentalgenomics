#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../log/1_sqanti.o
#SBATCH --error=../log/1_sqanti.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/SFARI_characterisation.config
source $SC_ROOT/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation

# directory
WK_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal
RB_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/Rerun


##-------------------------------------------------------------------------
mkdir -p ${WK_DIR}/1_SQANTI3Filtered
source activate sqanti2_py3

python $SQANTI3_DIR/sqanti3_filter.py rules ${RB_DIR}/test2_classification.txt --faa=${RB_DIR}/test2_corrected.faa --gtf=${RB_DIR}/test2_corrected.gtf \
-j=${relaxedJson} \
-o=merged \
-d=${WK_DIR}/1_SQANTI3Filtered \
--skip_report &> ${WK_DIR}/1_SQANTI3Filtered/merged_sqanti_filter.log

cd ${WK_DIR}/1_SQANTI3Filtered

# LOGEN: subset cupcake classification file by target genes
# merge cupcake classification file with abundance
subset_quantify_filter_tgenes.R \
--classfile ${WK_DIR}/1_SQANTI3Filtered/merged_RulesFilter_result_classification.txt \
--expression ${RB_DIR}/demux_fl_count.csv \
--target_genes ${TGENES_TXT} 

# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${WK_DIR}/1_SQANTI3Filtered/merged_RulesFilter_result_classification.txt \
--expression ${RB_DIR}/demux_fl_count.csv  \
--target_genes ${TGENES_TXT} \
--filter --nsample=2 --nreads=2

