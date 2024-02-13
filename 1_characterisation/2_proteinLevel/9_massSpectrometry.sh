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

#wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#wget https://zenodo.org/record/5076056/files/Task1SearchTaskconfig_orf.toml

refDir=/lustre/projects/Research_Project-MRC148213/lsl693/references
gencode_gtf=${refDir}/human/gencode.v40.annotation.gtf
gencode_fa=${refDir}/human/hg38.fa
gencode_transcript_fasta=${refDir}/human/gencode.v40.transcripts.fa
gencode_translation_fasta=${refDir}/human/gencode.v40.pc_translations.fa
uniprot_fasta=${refDir}/metamorpheus/uniprot_sprot.fasta

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
source activate mercury
export PATH=$PATH:/lustre/home/sl693/.dotnetexport PATH=$PATH:/lustre/home/sl693/.dotnet

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

#source /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteogenomics_functions.sh
#source /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteomics.config
#run_hybrid_annotation

Part1WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/8_longReadProteogenomics/longReadProteogenomics
Part2WKD=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/13_massSpec
run_metamorpheus $Part1WKD/7_classified_protein/$name".filtered_protein.fasta" $mass_spec $Part2WKD $name"_filtered"
#run_metamorpheus $protein_wkd/$name".protein_refined.fasta" $protein_dir $protein_wkd/All $name"_refined"
#run_metamorpheus $protein_wkd/$name"_hybrid.fasta" $protein_dir $protein_wkd/All $name"_hybrid"
#run_metamorpheus $WKD/gencode_protein.fasta $protein_dir $protein_wkd/All Gencode
#run_metamorpheus $uniprot_fasta $protein_dir $protein_wkd/All UniProt


