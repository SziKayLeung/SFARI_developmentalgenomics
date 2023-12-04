#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../log/1_sqanti.o
#SBATCH --error=../log/1_sqanti.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/1_transcriptLevel
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE
source $SC_ROOT/SFARI_characterisation.config
source $SC_ROOT/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${FICLE_ROOT}
export PATH=$PATH:${FICLE_ROOT}/reference

# directory
WK_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal
RB_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/FINAL
RBWhole_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/11_cupcake/Rerun

##-------------------------------------------------------------------------
# LOGEN: subset cupcake classification file by target genes
# merge cupcake classification file with abundance
source activate nanopore
subset_quantify_filter_tgenes.R \
--classfile ${WK_DIR}/1_SQANTI3/SQANTI3_collapse_options_RulesFilter_result_classification.txt \
--expression ${RB_DIR}/demux_fl_count.csv \
--target_genes ${TGENES_TXT} 

# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${WK_DIR}/1_SQANTI3/SQANTI3_collapse_options_RulesFilter_result_classification.txt \
--expression ${RB_DIR}/demux_fl_count.csv  \
--target_genes ${TGENES_TXT} \
--filter --nsample=2 --nreads=2

# for whole transcriptome dataset
#mkdir -p ${WK_DIR}/WholeTranscriptome
#cp ${RBWhole_DIR}/SQANTI3_whole_RulesFilter_result_classification.txt ${WK_DIR}/WholeTranscriptome
#subset_quantify_filter_tgenes.R \
#--classfile ${WK_DIR}/WholeTranscriptome/SQANTI3_whole_RulesFilter_result_classification.txt \
#--expression ${RBWhole_DIR}/demux_fl_count.csv \
#--filter --nsample=2 --nreads=2

# working variables
collapsedGtf=${WK_DIR}/1_SQANTI3/sqantiqc_collapsed_corrected.gtf
collapsedFasta=${WK_DIR}/1_SQANTI3/sqantiqc_collapsed_corrected.fasta
finalanno=${WK_DIR}/1_SQANTI3Filtered/test/SQANTI3_collapse_options_RulesFilter_result_classification.targetgenes_counts_filtered.txt
finaliso=${WK_DIR}/1_SQANTI3Filtered/test/SQANTI3_collapse_options_RulesFilter_result_classification.targetgenes_filtered_isoforms.txt
prefix=sqantiqc_collapsed_corrected_counts_filtered

# LOGEN: subset fasta and gtf using the finalised list of target gene isoforms
source activate sqanti2_py3
subset_fasta_gtf.py --gtf ${collapsedGtf} -i ${finaliso} -o counts_filtered -d ${WK_DIR}/1_SQANTI3
subset_fasta_gtf.py --fa ${collapsedFasta} -i ${finaliso} -o counts_filtered -d ${WK_DIR}/1_SQANTI3

gtfToGenePred ${prefix}.gtf ${prefix}.genePred
genePredToBed ${prefix}.genePred > ${prefix}.bed12
sort -k1,1 -k2,2n ${prefix}.bed12 > ${prefix}_sorted.bed12


##-------------------------------------------------------------------------

# run CPAT
mkdir -p ${WK_DIR}/2_cpat; cd ${WK_DIR}/2_cpat
cpat.py -x $HEXAMER -d $LOGITMODEL -g ${WK_DIR}/1_SQANTI3/${prefix}.fa --antisense --min-orf=50 --top-orf=50 -o ${prefix} 2> ${prefix}"_cpat.e"

# extract_best_orf <sample> <root_dir>
cpatRoot=${WK_DIR}/2_cpat/${prefix}
extract_fasta_bestorf.py --fa ${cpatRoot}".ORF_seqs.fa" --orf ${cpatRoot}".ORF_prob.best.tsv" --o_name merged_bestORF --o_dir ${WK_DIR}/2_cpat &> orfextract.log


##-------------------------------------------------------------------------

# run FICLE
mkdir -p ${WK_DIR}/3_ficle ${WK_DIR}/3_ficle/references
cd ${WK_DIR}/3_ficle
subset_reference_by_gene.py --r ${GENOME_GTF} --glist MECP2 --o ${WK_DIR}/3_ficle/references

inputRef=${WK_DIR}/3_ficle/references
inputGtf=${WK_DIR}/1_SQANTI3Filtered/merged.filtered_counts_filtered.gtf
inputBed=${WK_DIR}/1_SQANTI3Filtered/merged.filtered.sorted.bed12
inputClass=${WK_DIR}/1_SQANTI3Filtered/merged_RulesFilter_result_classification.targetgenes_counts_filtered.txt
inputOrf=${WK_DIR}/2_cpat/merged.ORF_prob.best.tsv
outputDir=${WK_DIR}/3_ficle

ficle.py --gene=MECP2 -r=${inputRef} -b=${inputBed} -g=${inputGtf} -c=${inputClass} --cpat=${inputOrf} -o=${outputDir} &> MECP2_characterise.log

