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
#SBATCH --output=9_massSpectrometry_3.o
#SBATCH --error=9_massSpectrometry_3.e

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

for SLURM_ARRAY_TASK_ID in {21..34}; do 
  echo $SLURM_ARRAY_TASK_ID
  
  
  protein=${protein_data[${SLURM_ARRAY_TASK_ID}]}
  echo $protein
  
  # run_metamorpheus <input_fasta>
  run_metamorpheus(){
    # variables 
    input_fasta=$1
    protein_dir=$2
    output_dir=$3
    output_name=$4
    
    echo "Processing $input_fasta"
    protein_data=$(for i in $protein_dir/*raw*; do echo $i; done)
    echo $protein_data
    
    cd $output_dir 
    dotnet $metamorpheus_dir/CMD.dll -g -o ./toml --mmsettings ./settings
    dotnet $metamorpheus_dir/CMD.dll -d $input_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $protein_data -t $metamorpheus_toml -v normal --mmsettings settings -o ./$output_name"_search_results"
    mv $output_name"_search_results"/Task1SearchTask/AllPeptides.psmtsv $output_name"_search_results"/Task1SearchTask/AllPeptides.$output_name".psmtsv"
    mv $output_name"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.tsv $output_name"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.$output_name".filtered.tsv"  
  }
  
  convert_gtf_bed12(){
    echo "Processing $1.gtf for conversion to bed12"
    gtfToGenePred $1.gtf $1.genePred
    genePredToBed $1.genePred $1.bed12
    if [ $3 == "make_region" ]; then
    # Squish intronic regions
    python $LREAD/visualization_track/src/make_region_bed_for_ucsc.py --name $2 --sample_gtf $1.gtf --reference_gtf $GENOME_GTF
    fi
  }
  
  
  SFARIfasta=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics/7_classified_protein/Whole.filtered_protein.fasta
  Part1WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics
  Part2WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/13_massSpec
  mkdir -p ${Part2WKD}/input_massspec 
  
  cd ${Part2WKD}
  mkdir -p ${SLURM_ARRAY_TASK_ID}
  
  cd ${SLURM_ARRAY_TASK_ID}
  #yes y | dotnet $metamorpheus_dir/CMD.dll -g -o ./toml --mmsettings ./settings
  #yes y | dotnet $metamorpheus_dir/CMD.dll -d $SFARIfasta settings/Contaminants/MetaMorpheusContaminants.xml -s ${protein} -t $metamorpheus_toml -v normal --mmsettings settings -o ./search_results
  
  #SFARIfasta=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics/7_classified_protein/Whole.filtered_protein.fasta
  #dotnet $metamorpheus_dir/CMD.dll -g -o ./toml --mmsettings ./settings
  #dotnet $metamorpheus_dir/CMD.dll -d $SFARIfasta settings/Contaminants/MetaMorpheusContaminants.xml -s ${Part2WKD}/input_massspec/20210108_AD_1.raw  -t $metamorpheus_toml -v normal --mmsettings settings -o ./"SFARI_search_results"
  
  
  export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
  export LREAD=${SOFTDIR}/Long-Read-Proteogenomics/modules/
  GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.gtf
  
  mkdir -p visualisation novel_peptides
  cd visualisation
  python $LREAD/visualization_track/src/make_peptide_gtf_file.py \
  --name Sfari \
  --sample_gtf ${Part1WKD}/7_classified_protein/Whole_with_cds_filtered.gtf \
  --reference_gtf $GENOME_GTF \
  --peptides ${Part2WKD}/$SLURM_ARRAY_TASK_ID/search_results/Task1SearchTask/AllPeptides.psmtsv \
  --pb_gene ${Part1WKD}/7_classified_protein/Whole_genes.tsv \
  --gene_isoname ${Part1WKD}/3_reference_tables/gene_isoname.tsv \
  --refined_fasta $Part1WKD/7_classified_protein/Whole.protein_refined.fasta
  
  gtfToGenePred Sfari_peptides.gtf Sfari_peptides.genePred
  genePredToBed Sfari_peptides.genePred > Sfari_peptides.bed12
  python $LREAD/visualization_track/src/make_region_bed_for_ucsc.py --name Sfari --sample_gtf Sfari_peptides.gtf --reference_gtf $GENOME_GTF
  
  cd ../novel_peptides
  python $LREAD/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
  --pacbio_peptides ${Part2WKD}/$SLURM_ARRAY_TASK_ID/search_results/Task1SearchTask/AllPeptides.psmtsv \
  --gencode_fasta $Part1WKD/3_reference_tables/gencode_protein.fasta \
  --uniprot_fasta $uniprot_fasta \
  --name SFARI
done