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
#SBATCH --output=2_wholedataset_DIU-%A_%a.o
#SBATCH --error=2_wholedataset_DIU-%A_%a.e
#SBATCH --array=0-1


# 13/10/2023: run DIU in targeted group dataset
# 02/11/2023: run DIU in all datasets


##-------------------------------------------------------------------------

# source function script

module load Miniconda2
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/differential_analysis

# directory paths
SQANTIDIR=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3
DESEQDIR=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/7_differential/2_DTE
DIUDIR=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/7_differential/3_DIU
mkdir -p $DIUDIR/targeted $DIUDIR/whole
mkdir -p $DIUDIR/targeted/sex $DIUDIR/targeted/group $DIUDIR/whole/sex $DIUDIR/whole/group $DIUDIR/whole/allGroup $DIUDIR/whole/allSex

# classification file and target genes
SQANTIClassFile=${SQANTIDIR}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic.txt
METADIR=/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/7_differential/3_DIU/metadata
  
WholeGroupPhenotype=${METADIR}/WholeGroupPhenotype.csv
WholeGroupFactor=${METADIR}/WholeGroupFactors.txt

WholeSexPhenotype=${METADIR}/WholeSexPhenotype.csv
WholeSexFactor=${METADIR}/WholeSexFactors.txt

ProteinCodingGenes=${METADIR}/WholeProteinCodingGenes.txt
lncRNAGenes=${METADIR}/WholelncRNAGenes.txt

GeneTypes=($ProteinCodingGenes $lncRNAGenes)
GeneType=${GeneTypes[${SLURM_ARRAY_TASK_ID}]}


##-------------------------------------------------------------------------

# run DIU through each gene list
source activate edgeR
while read gene; do
  echo "$gene"
  # Whole group
  python ${LOGEN_ROOT}/differential_analysis/subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_whole_group_removinglateprenatal_TEST.csv -d ${DIUDIR}/whole/allGroup -I ${gene} -p ${WholeGroupPhenotype} -f ${WholeGroupFactor} -n "Whole"
  # Whole sex
  python ${LOGEN_ROOT}/differential_analysis/subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/output_norm_whole_sex_removinglateprenatal_TEST.csv -d ${DIUDIR}/whole/allSex -I ${gene} -p ${WholeSexPhenotype} -f ${WholeSexFactor} -n "Whole"
done <${GeneType}