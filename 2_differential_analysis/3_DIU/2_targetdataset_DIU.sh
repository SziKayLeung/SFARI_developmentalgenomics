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
#SBATCH --output=2b_targetedataset_DIU.o
#SBATCH --error=2b_targetedataset_DIU.o


# 13/10/2023: run DIU in targeted group dataset
# 02/11/2023: run DIU in all datasets


##-------------------------------------------------------------------------

# source function script

module load Miniconda2
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/differential_analysis

# directory paths
SQANTIDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI
DESEQDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/4_deseq2
DIUDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu
mkdir -p $DIUDIR/targeted $DIUDIR/whole
mkdir -p $DIUDIR/targeted/sex $DIUDIR/targeted/group $DIUDIR/whole/sex $DIUDIR/whole/group

# classification file and target genes
SQANTIClassFile=${SQANTIDIR}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples_classification_noMonoIntergenic_counts.txt
TGENETXT=/gpfs/mrc0/projects/Research_Project-MRC148213/vc362/fetalBrain/genes.txt
METADIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/5_diu/metadata/
  
  TargetedGroupPhenotype=${METADIR}/TargetedDIUphenotype.csv
TargetedGroupFactor=${METADIR}/TargetedGroupFactors.txt

TargetedSexPhenotype=${METADIR}/TargetedSexPhenotype.csv
TargetedSexFactor=${METADIR}/TargetedSexFactors.txt

WholeGroupPhenotype=${METADIR}/WholeGroupPhenotype.csv
WholeGroupFactor=${METADIR}/WholeGroupFactors.txt

WholeSexPhenotype=${METADIR}/WholeSexPhenotype.csv
WholeSexFactor=${METADIR}/WholeSexFactors.txt


##-------------------------------------------------------------------------

# run DIU through each target gene
source activate edgeR
while read gene; do
echo "$gene"
# Targeted group
#subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_targeted_group_removinglateprenatal.csv -d ${DIUDIR}/targeted/group -I ${gene} -p ${TargetedGroupPhenotype} -f ${TargetedGroupFactor} -n "Targeted"
# Targeted sex
subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_targeted_sex_removinglateprenatal.csv -d ${DIUDIR}/targeted/sex -I ${gene} -p ${TargetedSexPhenotype} -f ${TargetedSexFactor} -n "Targeted"
# Whole group
subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_whole_group_removinglateprenatal_TEST.csv -d ${DIUDIR}/whole/group -I ${gene} -p ${WholeGroupPhenotype} -f ${WholeGroupFactor} -n "Whole"
# Whole sex
subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_whole_sex_removinglateprenatal_TEST.csv -d ${DIUDIR}/whole/sex -I ${gene} -p ${WholeSexPhenotype} -f ${WholeSexFactor} -n "Whole"
done <${TGENETXT}