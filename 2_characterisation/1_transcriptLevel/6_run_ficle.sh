
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/sl693/scripts/FICLE
export PATH=$PATH:${FICLE_ROOT}

SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/
iRef=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references/annotation/gencode.v40.annotation.gtf
iBed12=${SQANTI}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.bed12
iGtf=${SQANTI}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf
iClass=${SQANTI}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt
iCpat=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.ORF_prob.best.tsv
oDir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle/


ficle.py --gene=GRIA3 --reference=${iRef} --input_bed=${iBed12} --input_gtf=${iGtf} --input_class=${iClass} --cpat=${iCpat} --output_dir=${oDir}

grep ONTX_7115_1001 ${iGtf} | head
grep ONTX_7115_1001 ${iClass} | head
