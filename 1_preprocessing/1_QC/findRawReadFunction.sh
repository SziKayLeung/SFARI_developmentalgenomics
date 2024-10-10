find_raw_read ONT1_12_2049
find_raw_read(){
  
  echo "Finding the raw read ID"
  cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/2_trimmed/targeted/duplicated
  for i in *duplicated_read_stats_final.csv*; do
  if grep -q $1 "$i"; then
  ontReadID=$(grep $1 "$i" | awk '{print $1}')
  fi
  done
  
  echo "Finding the sample with that raw read ID"
  cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/2_trimmed/targeted/duplicated/samples
  for i in *; do
  if grep -q "$ontReadID" "$i"; then
  sampleMatched="${i%_duplicated_reads.csv}"
  fi
  done
  
  echo "Finding the sample uploaded to SRA"
  manifest=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/WholeTargetedphenotype_manifest.csv
  uploadedSample=$(grep $sampleMatched $manifest | awk -F "," '{ print $1 }')
  
  echo "Extracting the raw read"
  SRA=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined/
    zgrep -A 3 $ontReadID $SRA/$uploadedSample.fastq.gz
  
}