#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=120:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-19 # 20 files split across 17,000 protein coding genes
#SBATCH --output=2_wholedataset_proteincoding_DIU-%A_%a.o
#SBATCH --error=2_wholedataset_proteincoding_DIU-%A_%a.e


# 30/01/2024: perform DIU on protein-coding genes from whole dataset
# 06/02/2024: repeat but with filtered SQANTI3 classification after removal of mono-exonic transcripts within multiexonic genes
# 05/11/2024: repeat after re-running whole+targeted analyis

##-------------------------------------------------------------------------

# source function script

module load Miniconda2
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/differential_analysis

# directory paths
SQANTIDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final
DESEQDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/18_deseq/2_DTE
DIUDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/18_deseq/3_DIU
mkdir -p ${DIUDIR}
mkdir -p ${DIUDIR}/whole
mkdir -p $DIUDIR/whole/allGroup $DIUDIR/whole/allSex

# classification file and target genes
SQANTIClassFile=${SQANTIDIR}/sqantifiltered_monoexonicfiltered_2reads2samples_classification_finalversion.txt
METADIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/18_deseq/3_DIU/metadata
  
WholeGroupPhenotype=${METADIR}/WholeGroupPhenotype.csv
WholeGroupFactor=${METADIR}/WholeGroupFactors.txt

WholeSexPhenotype=${METADIR}/WholeSexPhenotype.csv
WholeSexFactor=${METADIR}/WholeSexFactors.txt

ProteinCodingGenes=${METADIR}/WholeProteinCodingGenes.txt
lncRNAGenes=${METADIR}/WholelncRNAGenes.txt

# remove empty lines 
#cd ${METADIR}/splitWholeGenes/
#for file in *; do sed -i '/^$/d' "$file"; done

# setting up SLURM
WholeGenesLists=($(ls ${METADIR}/splitWholeGenes/*))
WholeGenesList=${WholeGenesLists[${SLURM_ARRAY_TASK_ID}]}
#echo "${#WholeGenesList[@]}"

##-------------------------------------------------------------------------

# run DIU through each gene list
source activate edgeR
while read gene; do
  if [ -f "${DIUDIR}/whole/allGroup/${gene}_classification.txt" ]; then
    echo "$gene already processed for group"
  else
    echo "To process ${gene} for group"
    # Whole group
    python ${LOGEN_ROOT}/differential_analysis/subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/DESeq2_whole_development_normAll.csv -d ${DIUDIR}/whole/allGroup -I ${gene} -p ${WholeGroupPhenotype} -f ${WholeGroupFactor} -n "Whole"
  fi
  
  if [ -f "${DIUDIR}/whole/allSex/${gene}_classification.txt" ]; then
    echo "$gene already processed for sex"
  else
    echo "To process ${gene} for sex"
    # Whole sex
    python ${LOGEN_ROOT}/differential_analysis/subset_and_run_DIU.py ${SQANTIClassFile} ${DESEQDIR}/DESeq2_whole_sex_normAll.csv -d ${DIUDIR}/whole/allSex -I ${gene} -p ${WholeSexPhenotype} -f ${WholeSexFactor} -n "Whole"
  fi

done < ${WholeGenesList}
