#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=9b_rerunNovelPeptide.o
#SBATCH --error=9b_rerunNovelPeptide.e

#wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#wget https://zenodo.org/record/5076056/files/Task1SearchTaskconfig_orf.toml

refDir=/lustre/projects/Research_Project-MRC148213/lsl693/references
gencode_gtf=${refDir}/human/gencode.v40.annotation.gtf
gencode_fa=${refDir}/human/hg38.fa
gencode_transcript_fasta=${refDir}/human/gencode.v40.transcripts.fa
gencode_translation_fasta=${refDir}/human/gencode.v40.pc_translations.fa
uniprot_fasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/homo_sapiens_uniprot.fasta

hexamer=${refDir}/CPAT/Human_Hexamer.tsv
logit_model=${refDir}/CPAT/Human_logitModel.RData

sqantiDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti
classification=${sqantiDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Whole_2reads2samples_monomultirem.txt
sqanti_fasta=${sqantiDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_Whole_2reads2samples_nomonointergenic.fasta
sqanti_gtf=${sqantiDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_Whole_2reads2samples_nomonointergenic.gtf

metamorpheus_dir=/lustre/projects/Research_Project-MRC148213/lsl693/software/MetaMorpheus
metamorpheus_toml=${refDir}/metamorpheus/Task1SearchTaskconfig_orf.toml
meta_rescue_toml=${refDir}/metamorpheus/Task1SearchTaskconfig_rescue_resolve.toml

mass_spec=/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/1_raw/C_Proteomics
name=Whole
coding_score_cutoff=0.0
min_junctions_after_stop_codon=2
lower_kb=1
upper_kb=4
lower_cpm=3

module load Miniconda2 
source activate sqanti2_py3
export PATH=$PATH:/lustre/home/sl693/.dotnet
#protein_data=(${mass_spec}/20210108_AD_1.raw ${mass_spec}/20210108_AD_2.raw)
protein_data=(
  ${mass_spec}/20210108_AD_1.raw
  ${mass_spec}/20210108_AD_2.raw
  ${mass_spec}/20210112_AD_3.raw
  ${mass_spec}/20210112_AD_4.raw
  ${mass_spec}/20210115_AD_5.raw
  ${mass_spec}/20210115_AD_6.raw
  ${mass_spec}/20210120_AD_7.raw
  ${mass_spec}/20210330_AD_8.raw
  ${mass_spec}/20210328_AD_9.raw
  ${mass_spec}/20210325_AD_10.raw
  ${mass_spec}/20210324_AD_11.raw
  ${mass_spec}/20210308_AD_12.raw
  ${mass_spec}/20210108_AD_13.raw
  ${mass_spec}/20210108_AD_14.raw
  ${mass_spec}/20210112_AD_15.raw
  ${mass_spec}/20210112_AD_16.raw
  ${mass_spec}/20210115_AD_17.raw
  ${mass_spec}/20210115_AD_18.raw
  ${mass_spec}/20210120_AD_19.raw
  ${mass_spec}/20210330_AD_20.raw
  ${mass_spec}/20210328_AD_21.raw
  ${mass_spec}/20210325_AD_22.raw
  ${mass_spec}/20210324_AD_23.raw
  ${mass_spec}/20210308_AD_24.raw
  ${mass_spec}/20210110_AD_25.raw
  ${mass_spec}/20210108_AD_26.raw
  ${mass_spec}/20210112_AD_27.raw
  ${mass_spec}/20210112_AD_28.raw
  ${mass_spec}/20210115_AD_29.raw
  ${mass_spec}/20210115_AD_30.raw
  ${mass_spec}/20210120_AD_31.raw
  ${mass_spec}/20210330_AD_32.raw
  ${mass_spec}/20210328_AD_33.raw
  ${mass_spec}/20210325_AD_34.raw
  ${mass_spec}/20210324_AD_35.raw
)

for SLURM_ARRAY_TASK_ID in {0..34}; do 
  echo $SLURM_ARRAY_TASK_ID


  SFARIfasta=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics/7_classified_protein/Whole.filtered_protein.fasta
  SFARIgtf=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf
  Part1WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics
  Part2WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/13_massSpec
  export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
  export LREAD=${SOFTDIR}/Long-Read-Proteogenomics/modules/
  GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.gtf
  
  cd ${Part2WKD}/${SLURM_ARRAY_TASK_ID}/novel_peptides
  #python $LREAD/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
  #--pacbio_peptides ${Part2WKD}/$SLURM_ARRAY_TASK_ID/search_results/Task1SearchTask/AllPeptides.psmtsv \
  #--gencode_fasta $Part1WKD/3_reference_tables/gencode_protein.fasta \
  #--uniprot_fasta $uniprot_fasta \
  #--name SFARI
  awk '{ print $3 }' SFARI.pacbio_novel_peptides_to_uniprot.tsv > SFARI.pacbio_novel_peptides_seq.tsv
  grep -f SFARI.pacbio_novel_peptides_seq.tsv ../../AllSfari_peptide.bed12 > SFARI.pacbio_novel_peptides_seq.bed12
  sort SFARI.pacbio_novel_peptides_seq.bed12 | uniq > SFARI.pacbio_novel_peptides_seq_unique.bed12
  awk '{ print $2 }' SFARI.pacbio_novel_peptides_to_uniprot.tsv > SFARI.pacbio_novel_peptides_id.tsv
  grep -f SFARI.pacbio_novel_peptides_id.tsv ${SFARIgtf} >  SFARI.pacbio_novel_peptides.gtf

done