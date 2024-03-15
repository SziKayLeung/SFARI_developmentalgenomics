#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
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
module load R
#source activate edgeR

rootDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI
deseqRscript=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/2_differential_analysis/deseq_functions.R
inputPhenotype=${rootDir}/0_metadata/Wholephenotype_fixedsex.csv
inputExpression=${rootDir}/10_deseq/WholeTargeted_demux.csv
inputClassfile=${rootDir}/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples_monomultirem.txt


##-------------------------------------------------------------------------
# gene level
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="development" --level="gene" --dataset="whole" --directory=${rootDir}/10_deseq/1_DGE
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="sex" --level="gene" --dataset="whole" --directory=${rootDir}/10_deseq/1_DGE

# transcript level
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="development" --level="transcript" --dataset="whole" --directory=${rootDir}/10_deseq/2_DTE
Rscript ${deseqRscript} --phenotype=${inputPhenotype} --expression=${inputExpression} --classfile=${inputClassfile} --design="sex" --level="transcript" --dataset="whole" --directory=${rootDir}/10_deseq/2_DTE
