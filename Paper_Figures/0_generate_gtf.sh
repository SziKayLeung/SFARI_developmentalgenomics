
input=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final
output=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/specificGenes

grep -w ONT18.5258.1932 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT18.5258.1932.gtf
grep -w ONT10.5139.1910 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT10.5139.1910.gtf
grep -w ONT2.10213.11813 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT2.10213.11813.gtf