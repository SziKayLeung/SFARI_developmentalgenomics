#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --mem=200G
#SBATCH --output=1_demux4.o
#SBATCH --error=1_demux4.e

SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
TCLEANDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/BAMford
COLLAPSED=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed
SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
source $SC_ROOT/SFARI_characterisation.config

module load Miniconda2
source activate sqanti2_py3

#adapt_cupcake_to_ont.py ${TCLEANDIR} -o WholeTargeted -i clean.fa
adapt_cupcake_to_ont.py  ${TCLEANDIR} -s ${TCLEANDIR}/Whole1192133_aligned_clean.fa -o WholeTargeted -i clean.fa

cat ${COLLAPSED}/*WholeTargeted*read_stat.txt > ${COLLAPSED}/WholeTargeted_cleaned_aligned_merged_collapsed_allchrs.read_stat.txt
demux_cupcake_collapse.py ${COLLAPSED}/WholeTargeted_cleaned_aligned_merged_collapsed_allchrs.read_stat.txt ${TCLEANDIR}/WholeTargeted_sample_id.csv 

#source activate nanopore
#subset_quantify_filter_tgenes.R \
#--classfile ${SQANTI}/Targeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification.txt \
#--expression ${COLLAPSED}/demux_fl_count.csv \
#--target_genes ${TGENES_TXT} 