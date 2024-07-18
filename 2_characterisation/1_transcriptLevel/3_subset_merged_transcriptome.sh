#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
##SBATCH --output=../log/3_subset_merged_transcriptome.o
##SBATCH --error=../log/3_subset_merged_transcriptome.e
#SBATCH --mem=200G # specify bytes memory to reserve

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/1_transcriptLevel
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 

#wholeTargetedDir=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/
#wholeTargetedGff=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/merged_sorted.gff
#WK_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/

source activate sqanti2_py3
#Rscript ${SC_ROOT}/3_subset_merged_transcriptome.R
#subset_fasta_gtf.py ${wholeTargetedGff} --gtf -i ${wholeTargetedDir}/filteredIds.csv -o merged_sorted_filtered2reads2samples -d ${WK_DIR}

# 2 reads 2 samples, with no monointergenic reads
wholeTargetedDir=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI
wholeTargetedGff=${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.gtf
wholeTargetedFasta=${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.fasta
RetainedIDs=${wholeTargetedDir}/WholeTargeted_RulesFilter_2reads2samples_nomonointergenic_isoforms.txt 
subset_fasta_gtf.py ${wholeTargetedGff} --gtf -i ${RetainedIDs} -o 2reads2samples_nomonointergenic -d ${wholeTargetedDir}
subset_fasta_gtf.py ${wholeTargetedFasta} --fa -i ${RetainedIDs} -o 2reads2samples_nomonointergenic -d ${wholeTargetedDir}
seqtk subseq ${wholeTargetedFasta} ${RetainedIDs} > ${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_nomonointergenic.fasta

# whole dataset only 
wholeTargetedDir=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI
wholeTargetedGff=${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.gtf
wholeTargetedFasta=${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.fasta
RetainedIDs=${wholeTargetedDir}/WholeTargeted_RulesFilter_Whole_2reads2samples_nomonointergenic_isoforms.txt
subset_fasta_gtf.py ${wholeTargetedGff} --gtf -i ${RetainedIDs} -o Whole_2reads2samples_nomonointergenic -d ${wholeTargetedDir}
seqtk subseq ${wholeTargetedFasta} ${RetainedIDs} > ${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_Whole_2reads2samples_nomonointergenic.fasta


#grep -f ${RetainedIDs} ${wholeTargetedDir}/WholeTargeted_demux_2reads2samples_SQANTIfiltered.csv > ${wholeTargetedDir}/WholeTargeted_demux_2reads2samples_nomonointergenic_SQANTIfiltered.csv

#gtfToGenePred WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.genepred
#genePredToBed WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.genepred > WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.bed12
colour_transcripts_by_countandpotential.py \
 --a=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_demux_2reads2samples_SQANTIfiltered_dataset.csv \
 -s=human \
 --bed=${wholeTargetedDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.bed12 \
 --cpat=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.ORF_prob.best.tsv \
 --noORF=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.no_ORF.txt \
 --d=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/Tracks/ \
 --output=WholeTargeted 
 
 #paste0(dirnames$wholetarg_SQ, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples.txt")
 #subset_and_run_DIU.py WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_Targeted_2reads2samples_classification_noMonoIntergenic_counts.txt output_norm_targeted_group_removinglateprenatal.csv -I SEPTIN4
