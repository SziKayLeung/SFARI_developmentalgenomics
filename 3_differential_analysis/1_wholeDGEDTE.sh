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
#SBATCH --output=1_wholeDGEDTE.o
#SBATCH --error=1_wholeDGEDTE.e

## print start date and time
echo Job started on:
date -u

##-------------------------------------------------------------------------
module load Miniconda2
module load R/4.2.2-foss-2022b

#source activate edgeR

rootDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI
deseqRscript=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/3_differential_analysis/deseq_functions.R
inputPhenotype=${rootDir}/0_metadata/Wholephenotype_fixedsex.csv
inputExpression=${rootDir}/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron_classification.txt
inputClassfile=${rootDir}/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron_classification.txt
 
 
mkdir -p ${rootDir}/C_Whole_Targeted/18_deseq/1_DGE ${rootDir}/C_Whole_Targeted/18_deseq/2_DTE

##-------------------------------------------------------------------------
# gene level
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --classfile=${inputClassfile} --design="development" --level="gene" --dataset="whole" --directory=${rootDir}/C_Whole_Targeted/18_deseq/1_DGE
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --classfile=${inputClassfile} --design="sex" --level="gene" --dataset="whole" --directory=${rootDir}/C_Whole_Targeted/18_deseq/1_DGE

# transcript level
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="development" --level="transcript" --dataset="whole" --directory=${rootDir}/C_Whole_Targeted/18_deseq/2_DTE
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="sex" --level="transcript" --dataset="whole" --directory=${rootDir}/C_Whole_Targeted/18_deseq/2_DTE
