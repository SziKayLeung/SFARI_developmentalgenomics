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
#SBATCH --output=7_gffcompareDatasets.o
#SBATCH --error=7_gffcompareDatasets.e

# 30/10/2024: overlap between PacBio whole cortex dataset 
# 04/11/2024: rerun with gtf filtered by intergenicGenicIntron

module load Miniconda2/4.3.21

sqantiDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final
sqanti_whole_gtf=sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf
sqanti_wholetargeted_gtf=sqantifiltered_monoexonicfiltered_2reads2samples_intergenicGenicIntron.filtered.gtf

outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/14_OverlapDatasets
mkdir -p ${outputDir}/cellReports2021 ${outputDir}/directRNA ${outputDir}/BDRNatureComms2024 ${outputDir}/Patowary2024

## ----- run gffcompare on PacBio dataset ---- 

PacBioDir=/lustre/projects/Research_Project-MRC148213/lsl693/PacBioPaper/SQANTI2/HumanCTX
PacBio_gtf=HumanCTX.collapsed_classification.filtered_lite.gtf

source activate sqanti2_py3 
cd ${outputDir}/cellReports2021
cp ${sqantiDir}/${sqanti_whole_gtf} .
cp ${PacBioDir}/${PacBio_gtf} .
PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r ${sqanti_whole_gtf} ${PacBio_gtf} -o sfari_PacBio
gffcompare -r ${PacBio_gtf} ${sqanti_whole_gtf} -o PacBio_sfari


## ----- run gffcompare on direct RNA dataset dataset ---- 

directRNADir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/dRNA/Rosie/9_sqanti_final
directRNA_gtf=sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf

source activate sqanti2_py3 
cd ${outputDir}/directRNA
cp ${sqantiDir}/${sqanti_whole_gtf} .
cp ${directRNADir}/${directRNA_gtf} .
PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r ${sqanti_whole_gtf} ${directRNA_gtf} -o sfari_dRNA
gffcompare -r ${directRNA_gtf} ${sqanti_whole_gtf} -o dRNA_sfari


## ----- run gffcompare on BDR targeted dataset ---- 

ADBDRDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3
BDR_gtf=ontBDR_collapsed.filtered_counts_filtered.gtf

source activate sqanti2_py3 
cd ${outputDir}/BDRNatureComms2024
cp ${sqantiDir}/${sqanti_wholetargeted_gtf} .
cp ${ADBDRDir}/${BDR_gtf} .
PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r ${sqanti_wholetargeted_gtf} ${BDR_gtf} -o sfari_BDR
gffcompare -r ${BDR_gtf} ${sqanti_wholetargeted_gtf} -o BDR_sfari


## ----- run gffcompare on Patowary 2024 dataset ---- 

source activate sqanti2_py3 
cd ${outputDir}/Patowary2024
cp ${sqantiDir}/${sqanti_wholetargeted_gtf} .
# uploaded Patowary dataset: https://github.com/gandallab/Dev_Brain_IsoSeq/blob/main/data/sqanti/cp_vz_0.75_min_7_recovery_talon_corrected.gtf.gz
gunzip cp_vz_0.75_min_7_recovery_talon_corrected.gtf.gz
# https://github.com/gandallab/Dev_Brain_IsoSeq/blob/main/data/sqanti/cp_vz_0.75_min_7_recovery_talon_classification.txt.gz
gunzip cp_vz_0.75_min_7_recovery_talon_classification.txt.gz
Patowary_gtf=cp_vz_0.75_min_7_recovery_talon_corrected.gtf

PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r ${sqanti_wholetargeted_gtf} ${Patowary_gtf} -o sfari_Patowary
gffcompare -r ${Patowary_gtf} ${sqanti_wholetargeted_gtf} -o Patowary_sfari