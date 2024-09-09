#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=log/merged_post_sqanti.o
#SBATCH --error=log/merged_post_sqanti.e


# 09/09/2024: post sqanti per chromosome

##-------------------------------------------------------------------------

# paths and variables
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/6_demux
SQNATI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/
MERGED_SQANTI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti_merged
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/target_gene_annotation

##-------------------------------------------------------------------------

cd ${MERGED_SQANTI_DIR}

# merge classification file
ls ${SQNATI_DIR}/*RulesFilter_result_classification.txt > WholeTargeted_RulesFilter_result_classification_filenames.txt

# FNR==1 && NR!=1 { next }: This skips the first row (FNR==1) for all files except the first (NR!=1) file
awk 'FNR==1 && NR!=1 { next } { print }' $(cat WholeTargeted_RulesFilter_result_classification_filenames.txt) > WholeTargeted_RulesFilter_result_classification.txt
awk '$NF == "Isoform"' WholeTargeted_RulesFilter_result_classification.txt > WholeTargeted_RulesFilter_result_classification_isoform.txt

# merge fasta
ls ${SQNATI_DIR}/*_corrected.fasta > WholeTargeted_correctedfasta_filenames.txt
cat $(cat WholeTargeted_correctedfasta_filenames.txt) > WholeTargeted_corrected.fasta

# merge gtf
ls ${SQNATI_DIR}/*corrected.gtf > WholeTargeted_correctedgtf_filenames.txt
cat $(cat WholeTargeted_correctedgtf_filenames.txt) > WholeTargeted_corrected.gtf

# merge gff
ls ${SQNATI_DIR}/*corrected.gtf.cds.gff > WholeTargeted_correctedgff_filenames.txt
cat $(cat WholeTargeted_correctedgff_filenames.txt) > WholeTargeted_corrected.gtf.cds.gff

# demux files
ls ${DEMUX_DIR}/*_fl_count.csv > WholeTargeted_fl_count_filenames.txt
awk 'FNR==1 && NR!=1 { next } { print }' $(cat WholeTargeted_fl_count_filenames.txt) > WholeTargeted_fl_count.csv

# filter by 2 reads, 2 samples
module load R
subset_quantify_filter_tgenes.R \
--classfile ${MERGED_SQANTI_DIR}/WholeTargeted_RulesFilter_result_classification_isoform.txt \
--expression ${MERGED_SQANTI_DIR}/WholeTargeted_fl_count.csv --nsample=2 --nreads=2

# filter in the whole and targeted dataset, transcripts that are in the off-reads but not seen in the whole dataset

