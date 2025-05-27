
input=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final
output=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final/specificGenes

grep -w ONT17.1148 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf > ${output}/CACNA1G.gtf
grep -w ONT7.4758 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf > ${output}/FOXP2.gtf
grep -w ONT11.3444 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf > ${output}/DAGLA.gtf
grep -w ONT2.3331 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples.filtered.gtf > ${output}/RPS27A.gtf
grep -w ONT14.1780 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/DLGAP5.gtf
grep -w ONT4.13313 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/GPM6A.gtf
grep -w ONTX.6284 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/MORF4L2.gtf
grep -w ONT20.3125 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/GNAS.gtf
grep -w ONT18.5258.1932 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT18.5258.1932.gtf
grep -w ONT10.5139.1910 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT10.5139.1910.gtf
grep -w ONT2.10213.11813 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT2.10213.11813.gtf
grep -w ONT12.2697.32229 ${input}/sqantifiltered_monoexonicfiltered_2reads2samples_whole_intergenicGenicIntron.filtered.gtf > ${output}/ONT12.2697.32229.gtf