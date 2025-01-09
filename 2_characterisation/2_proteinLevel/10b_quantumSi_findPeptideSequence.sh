#!/bin/bash

# find the peptide sequence on the top 2 and 20 SFARI peptides predicted to be differentially expressed between prenatal and postnatal
# IDs from 10a_quantumSi_identify_DTP.R

cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/10_longReadProteogenomics/7_classified_protein
grep ONT18.5258.10189 Whole.protein_refined.fasta -A 10
grep ONT18.5258.10189 Whole.protein_refined.fasta -A 10

TopRankedNovelProteins=(
"ONTX.9761.5829"  
"ONT8.85.5521"    
"ONT5.1507.7410"  
"ONT15.230.3905"  
"ONTX.9761.6751"  
"ONT18.1124.5586" 
"ONT16.513.29803" 
"ONT16.73.106515"
"ONTX.9761.5886"  
"ONT4.10212.2149" 
"ONT4.10212.2149" 
"ONT2.9952.3158"  
"ONT4.10212.1410"
"ONT2.3644.280"   
"ONT8.2327.10140" 
"ONT2.3644.568"  
"ONT1.3388.6453"  
"ONTX.9761.6761" 
"ONT9.342.21595"  
"ONT3.117.12367"  
"ONTX.3093.14908")

for i in ${TopRankedNovelProteins[@]}; do 
  echo $i
  grep $i Whole.protein_refined.fasta -A 10
done