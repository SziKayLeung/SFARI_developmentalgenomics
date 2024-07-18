module load Miniconda2/4.3.21

source activate sqanti2_py3
SQANTI3_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/SQANTI3
LOGEN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
RB_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/Rerun
defaultJson=$SQANTI3_DIR/utilities/filter/filter_default.json
reducedJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
relaxedJson=$SQANTI3_DIR/utilities/filter/filter_relaxed.json

cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal
#cp ${RB_DIR}/test3_classification.txt . # usage of intropolis data
#cp /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0059_20220813_10780/Batch1/20220813_1259_3G_PAM33351_84e820b3/cupcake/Rerun/EMX1* .
#cp $RB_DIR/demux_fl_count.csv .


sed '1p;/ENSG00000135638.14/!d' $RB_DIR/merged_mapped_classification.txt > EMX1_classification.txt
awk -F '\t' '{print $1}' EMX1_classification.txt | sed '1d' > EMX1_isoform_id.txt
(head -n 1 demux_fl_count.csv && grep -F -f EMX1_isoform_id.txt demux_fl_count.csv) > EMX1_demux_fl_count.csv
(head -n 1 $RB_DIR/merged_mapped_junctions.txt && grep -F -f EMX1_isoform_id.txt $RB_DIR/merged_mapped_junctions.txt) > EMX1_junctions.txt

python $LOGEN/miscellaneous/subset_fasta_gtf.py --gtf $RB_DIR/merged_mapped_corrected.gtf -i EMX1_isoform_id.txt -d . -o EMX1
python $LOGEN/miscellaneous/subset_fasta_gtf.py --fa $RB_DIR/merged_mapped_corrected.fasta -i EMX1_isoform_id.txt -d . -o EMX1

python $SQANTI3_DIR/sqanti3_filter.py rules EMX1_classification.txt \
--faa=merged_mapped_corrected_EMX1.fa \
--gtf=merged_mapped_corrected_EMX1.gtf \
-j=${reducedJson} \
-o=EMX1ReducedCoverage \
--skip_report &> EMX1ReducedCoverage_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules EMX1_classification.txt \
--faa=merged_mapped_corrected_EMX1.fa \
--gtf=merged_mapped_corrected_EMX1.gtf \
-j=${relaxedJson} \
-o=EMX1Relaxed \
--skip_report &> EMX1Relaxed_sqanti_filter.log

python $LOGEN/miscellaneous/subset_fasta_gtf.py --gtf $RB_DIR/merged_mapped_corrected.gtf -d . -o EMX1Removed -I PB.20140.14,PB.20140.20,PB.20140.21,PB.20140.22,PB.20140.23,PB.20140.24,PB.20140.8 

# usage of intropolis data
#sed '1p;/EMX1/!d' $RB_DIR/SQANTI3_test6_classification.txt > EMX1_intropolisClassification.txt
#awk -F '\t' '{print $1}' EMX1_intropolisClassification.txt | sed '1d' > EMX1_intrpolisIsoform_id.txt
#python $LOGEN/miscellaneous/subset_fasta_gtf.py --gtf $RB_DIR/SQANTI3_test6_corrected.gtf-i EMX1_intropolisIsoform_id.txt -d . -o EMX1_intropolis
#python $LOGEN/miscellaneous/subset_fasta_gtf.py --fa $RB_DIR/SQANTI3_test6_corrected.fasta -i EMX1_intropolisIsoform_id.txt -d . -o EMX1_intropolis

python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test6_classification.txt \
--faa=SQANTI3_test6_corrected.fasta \
--gtf=SQANTI3_test6_corrected.gtf \
-j=${defaultJson} \
-o=EMX1DefaultIntropolis \
--skip_report &> EMX1DefaultIntropolis_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test6_classification.txt \
--faa=SQANTI3_test6_corrected.fasta \
--gtf=SQANTI3_test6_corrected.gtf \
-j=$reducedJson \
-o=EMX1ReducedIntropolis \
--skip_report &> EMX1ReducedIntropolis_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test6_classification.txt \
--faa=SQANTI3_test6_corrected.fasta \
--gtf=SQANTI3_test6_corrected.gtf \
-j=${reducedJson} \
-o=EMX1RelaxedIntropolis \
--skip_report &> EMX1RelaxedIntropolis_sqanti_filter.log

# default parameters
python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test7_classification.txt \
--faa=SQANTI3_test7_corrected.fasta \
--gtf=SQANTI3_test7_corrected.gtf \
-j=${defaultJson} \
-o=EMX1Default \
--skip_report &> EMX1Default_sqanti_filter.log

#python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test7_classification.txt \
#--faa=SQANTI3_test7_corrected.fasta \
#--gtf=SQANTI3_test7_corrected.gtf \
#-j=$reducedJson \
#-o=EMX1Reduced \
#--skip_report &> EMX1Reduced_sqanti_filter.log

#python $SQANTI3_DIR/sqanti3_filter.py rules SQANTI3_test7_classification.txt \
#--faa=SQANTI3_test7_corrected.fasta \
#--gtf=SQANTI3_test7_corrected.gtf \
#-j=${relaxedJson} \
#-o=EMX1Relaxed \
#--skip_report &> EMX1Relaxed_sqanti_filter.log


python $SQANTI3_DIR/sqanti3_filter.py rules DDX5Intropolis_classification.txt \
--faa=DDX5Intropolis_corrected.fasta \
--gtf=DDX5Intropolis_corrected.gtf \
-j=${defaultJson} \
-o=DDX5Intropolis \
--skip_report &> DDX5Intropolis_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules DDX5NoIntropolis_classification.txt \
--faa=DDX5NoIntropolis_corrected.fasta \
--gtf=DDX5NoIntropolis_corrected.gtf \
-j=${defaultJson} \
-o=DDX5NoIntropolis \
--skip_report &> DDX5NoIntropolis_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules DDX5Intropolis_classification.txt \
--faa=DDX5Intropolis_corrected.fasta \
--gtf=DDX5Intropolis_corrected.gtf \
-j=${relaxedJson} \
-o=DDX5RelaxedIntropolis \
--skip_report &> DDX5RelaxedIntropolis_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules DDX5NoIntropolis_classification.txt \
--faa=DDX5NoIntropolis_corrected.fasta \
--gtf=DDX5NoIntropolis_corrected.gtf \
-j=${relaxedJson} \
-o=DDX5RelaxedNoIntropolis \
--skip_report &> DDX5RelaxedNoIntropolis_sqanti_filter.log

#### DDX5

sed '1p;/ENSG00000108654.16/!d' $RB_DIR/merged_mapped_classification.txt > DDX5_classification.txt
awk -F '\t' '{print $1}' DDX5_classification.txt | sed '1d' > DDX5_isoform_id.txt
(head -n 1 demux_fl_count.csv && grep -F -f DDX5_isoform_id.txt demux_fl_count.csv) > DDX5_demux_fl_count.csv
(head -n 1 $RB_DIR/merged_mapped_junctions.txt && grep -F -f DDX5_isoform_id.txt $RB_DIR/merged_mapped_junctions.txt) > DDX5_junctions.txt

python $LOGEN/miscellaneous/subset_fasta_gtf.py --gtf $RB_DIR/merged_mapped_corrected.gtf -i DDX5_isoform_id.txt -d . -o DDX5
python $LOGEN/miscellaneous/subset_fasta_gtf.py --fa $RB_DIR/merged_mapped_corrected.fasta -i DDX5_isoform_id.txt -d . -o DDX5

python $SQANTI3_DIR/sqanti3_filter.py rules DDX5_classification.txt \
--faa=merged_mapped_corrected_DDX5.fa \
--gtf=merged_mapped_corrected_DDX5.gtf \
-j=${reducedJson} \
-o=DDX5ReducedCoverage \
--skip_report &> DDX5ReducedCoverage_sqanti_filter.log

python $SQANTI3_DIR/sqanti3_filter.py rules DDX5_classification.txt \
--faa=merged_mapped_corrected_DDX5.fa \
--gtf=merged_mapped_corrected_DDX5.gtf \
-j=${relaxedJson} \
-o=DDX5Relaxed \
--skip_report &> DDX5Relaxed_sqanti_filter.log

python $LOGEN/miscellaneous/subset_fasta_gtf.py --gtf $RB_DIR/merged_mapped_corrected.gtf -d . -o DDX5Removed -i DDX5_noncanon_ID.txt
