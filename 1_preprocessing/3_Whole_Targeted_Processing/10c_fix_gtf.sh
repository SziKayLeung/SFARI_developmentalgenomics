#!/bin/bash

for i in {1..22} X Y; 
do 
	grep -w chr${i} /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_corrected.gtf.cds.gff > chr${i}.gff; 
	sed -i "s/novelGene_/novelGene_chr${i}_/g" chr${i}.gff;
#	cat chr${i}.gff >> /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_corrected_finalversion.gtf.cds.gff
done

Rscript /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/filter_gtf.R

for i in {1..22} X Y;
do
	tail -n+4 chr${i}.gff >> /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/sqantifiltered_monoexonicfiltered_2reads2samples_corrected_finalversion.gtf.cds.gff
done
