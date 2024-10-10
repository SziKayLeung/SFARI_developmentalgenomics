#!/bin/bash

# paths
pbmm2Aligned=/media/disk2/sfari_RB/5_cupcake/5_align/
DEMUX_DIR=/media/disk2/sfari_RB/7_demux/
export PATH=$PATH:~/installs/LOGen/assist_ont_processing/

mamba activate SQANTI3.env

##-------------------------------------------------------------------------

# generate sample_id.csv for targeted and whole dataset
adapt_cupcake_to_ont.py ${pbmm2Aligned} -o Targeted -i mapped_filtered_sorted.fasta -d ${DEMUX_DIR}
adapt_cupcake_to_ont.py ${pbmm2Aligned} -o Whole -i aligned_clean_aligned_filtered.fasta -d ${DEMUX_DIR}

# remoe "combined_sorted" string in Targeted_sample_id.csv (due to the input of the fasta files)
cd ${DEMUX_DIR}
awk -F',' '{gsub("_combined_sorted", "", $2); print $1 "," $2}' Targeted_sample_id.csv > Targeted_sample_id_mod.csv

# concatenate whole and targeted sample_id.csv for downstream demux
cat Whole_sample_id.csv Targeted_sample_id_mod.csv  > WholeTargeted_sample_id.csv
