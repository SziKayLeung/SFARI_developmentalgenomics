# Differential expression & Splicing analysis 

Analyses were performed on **tappAS** (java) to identify differentially expressed and spliced transcripts.

**Input files**:
1. SQANTI3 gff file 
2. Expression matrix from TALON (ONT FL reads)
3. Design file 

Plots of normalised expression were generated using [custom scripts](https://github.com/SziKayLeung/SFARI_developmentalgenomics/blob/master/2_differential_analysis/0_source_differential_functions.R) and output normalised expression matrix from tappAS (stored in [config.file](https://github.com/SziKayLeung/SFARI_developmentalgenomics/blob/master/2_differential_analysis/sfari_differential.config)).  

## Analysis
1. Adult vs Fetal: [Case-control analysis](https://github.com/SziKayLeung/SFARI_developmentalgenomics/blob/master/2_differential_analysis/2_post_tappas_age.R)
2. Adult vs Fetal, Male vs Female: [Multiple Series Time-Course analysis](https://github.com/SziKayLeung/SFARI_developmentalgenomics/blob/master/2_differential_analysis/3_post_tappas_agevssex.R) 