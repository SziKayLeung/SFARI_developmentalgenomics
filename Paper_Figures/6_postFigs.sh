
outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/0_output
sqantiDir=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/

cd $outputDir
grep -f uniqueWholeIso.txt ${sqantiDir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf > uniqueWholedataset.gtf
