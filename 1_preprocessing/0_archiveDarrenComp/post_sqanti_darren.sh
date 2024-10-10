#!/bin/bash

##-------------------------------------------------------------------------

# paths
UTILS_DIR=/home/darren/installs/SFARI_developmentalgenomics/1_preprocessing
MERGED_CHROM_DIR=/media/disk2/sfari_RB/5_cupcake/mergedChrom
SPLIT=/media/disk2/sfari_RB/6_sqanti/split
SQANTI_OUTPUT=/media/disk2/sfari_RB/6_sqanti/sqanti_relaxed
MERGED_SQANTI_DIR=/media/disk2/sfari_RB/6_sqanti/sqanti_relaxed/sqanti_merged

##-------------------------------------------------------------------------

mkdir -p ${MERGED_SQANTI_DIR}
cd ${MERGED_SQANTI_DIR}

# merge classification file
ls ${SQANTI_OUTPUT}/*RulesFilter_result_classification.txt > WholeTargeted_RulesFilter_result_classification_filenames.txt

# FNR==1 && NR!=1 { next }: This skips the first row (FNR==1) for all files except the first (NR!=1) file
awk 'FNR==1 && NR!=1 { next } { print }' $(cat WholeTargeted_RulesFilter_result_classification_filenames.txt) > WholeTargeted_RulesFilter_result_classification.txt
awk '$NF == "Isoform"' WholeTargeted_RulesFilter_result_classification.txt > WholeTargeted_RulesFilter_result_classification_isoform.txt

# merge fasta
ls ${SQANTI_OUTPUT}/*_corrected.fasta > WholeTargeted_correctedfasta_filenames.txt
cat $(cat WholeTargeted_correctedfasta_filenames.txt) > WholeTargeted_corrected.fasta

# merge gtf
ls ${SQANTI_OUTPUT}/*corrected.gtf > WholeTargeted_correctedgtf_filenames.txt
cat $(cat WholeTargeted_correctedgtf_filenames.txt) > WholeTargeted_corrected.gtf

# merge gff
ls ${SQANTI_OUTPUT}/*corrected.gtf.cds.gff > WholeTargeted_correctedgff_filenames.txt
cat $(cat WholeTargeted_correctedgff_filenames.txt) > WholeTargeted_corrected.gtf.cds.gff

# demux files
#ls ${DEMUX_DIR}/*_fl_count.csv > WholeTargeted_fl_count_filenames.txt
#awk 'FNR==1 && NR!=1 { next } { print }' $(cat WholeTargeted_fl_count_filenames.txt) > WholeTargeted_fl_count.csv

# filter by 2 reads, 2 samples
#module load R
#subset_quantify_filter_tgenes.R \
#--classfile ${MERGED_SQANTI_DIR}/WholeTargeted_RulesFilter_result_classification_isoform.txt \
#--expression ${MERGED_SQANTI_DIR}/WholeTargeted_fl_count.csv --nsample=2 --nreads=2

# filter in the whole and targeted dataset, transcripts that are in the off-reads but not seen in the whole dataset