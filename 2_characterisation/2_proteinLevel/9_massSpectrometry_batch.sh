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
#SBATCH --array 0-34 # 35 samples
#SBATCH --output=9_massSpectrometry-%A_%a.o
#SBATCH --error=9_massSpectrometry-%A_%a.e

# 7/11/2024: batch run alignment of mass-spectrometry 30 samples

#-----------------------------------------------------------------------

module load Miniconda2 
source activate sqanti2_py3
export PATH=$PATH:/lustre/home/sl693/.dotnet
source /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/2_characterisation/2_proteinLevel/0_mass_spec.config

sampleNum=$((SLURM_ARRAY_TASK_ID + 1))
mkdir -p ${Part2WKD}/${sampleNum} 
mkdir -p ${Part2WKD}/${sampleNum}/visualisation 
mkdir -p ${Part2WKD}/${sampleNum}/novel_peptides

#-----------------------------------------------------------------------

mass_spec_data=${mass_spec_dataset[${SLURM_ARRAY_TASK_ID}]}

echo "############## run metamorpheus"
echo "with mass_spec: ${mass_spec_data}"
cd ${Part2WKD}/${sampleNum}
yes y | dotnet $metamorpheus_dir/CMD.dll -g -o ./toml --mmsettings ./settings > ./settings.log 2>&1
yes y | dotnet $metamorpheus_dir/CMD.dll -d ${protein_fasta} settings/Contaminants/MetaMorpheusContaminants.xml -s ${mass_spec_data} -t $metamorpheus_toml -v normal --mmsettings settings -o ./search_results > ./metamorpheus.log 2>&1

echo "############## visualise peptides as gtf"
cd visualisation
# using output from proteogenomics pipeline
python $LREAD/visualization_track/src/make_peptide_gtf_file.py \
--name Sfari \
--sample_gtf ${Part1WKD}/7_classified_protein/Whole_with_cds_filtered.gtf \
--reference_gtf ${gencode_gtf} \
--peptides ${Part2WKD}/$sampleNum/search_results/Task1SearchTask/AllPeptides.psmtsv \
--pb_gene ${Part1WKD}/7_classified_protein/Whole_genes.tsv \
--gene_isoname ${Part1WKD}/3_reference_tables/gene_isoname.tsv \
--refined_fasta $Part1WKD/7_classified_protein/Whole.protein_refined.fasta

gtfToGenePred Sfari_peptides.gtf Sfari_peptides.genePred
genePredToBed Sfari_peptides.genePred > Sfari_peptides.bed12
python $LREAD/visualization_track/src/make_region_bed_for_ucsc.py --name Sfari --sample_gtf Sfari_peptides.gtf --reference_gtf ${gencode_gtf}

cd ../novel_peptides
python $LREAD/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
--pacbio_peptides ${Part2WKD}/$sampleNum/search_results/Task1SearchTask/AllPeptides.psmtsv \
--gencode_fasta $Part1WKD/3_reference_tables/gencode_protein.fasta \
--uniprot_fasta $uniprot_fasta \
--name SFARI