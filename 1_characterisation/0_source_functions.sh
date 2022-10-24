
convert_gtf_bed12(){
  
  # variables 
  output_dir="$(dirname $1)" 
  sample=${1%.gtf} # removes .gtf
  
  source activate sqanti2_py3
  cd $output_dir
  
  gtfToGenePred $1 $sample.genePred
  genePredToBed $sample.genePred > $sample.bed12
  sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
  rm $sample.genePred $sample.bed12
  # Rscript script.R <input.classfile> <input_output_dir of bed file> <prefix>
  #Rscript $GENERAL/annotate_uscs_tracks.R $SQANTI $WKD $sample
  #bedToBigBed -extraIndex=name $sample $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
  
  #source deactivate
  #./bedToBigBed -as=bedExample2.as -type=bed9+3 -extraIndex=name $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
}

#run_gff_compare <output_name> <ref_gtf> <2nd_gtf> <root_dir>
run_gff_compare(){
  
  source activate nanopore
  
  mkdir -p $4/1_gffcompare; cd $4/1_gffcompare
  cp $2 . 
  cp $3 . 
  
  ref1="$(basename -- $2)"
  ref2="$(basename -- $3)"
  
  source activate sqanti2_py3 
  gffcompare -r $ref1 $ref2 -o $1
  
}


# run_sqanti3 <sample> <gtf> <root_dir>
run_sqanti3(){
  
  source activate sqanti2_py3
  
  # variable 
  sample=$1
  gtf=$2
  
  mkdir -p $3/2_sqanti3; cd $3/2_sqanti3
  
  # sqanti qc
  echo "Processing Sample $sample for SQANTI3 QC"
  python $SQANTI3_DIR/sqanti3_qc.py -v
  echo $GENOME_GTF
  echo $GENOME_FASTA
 
  echo "Processing basic commands"
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report pdf &> $sample.sqanti.qc.log
  
  echo "Processing Sample $sample for SQANTI2 filter"
  python $SQANTI3_DIR/sqanti3_RulesFilter.py $sample"_classification.txt" $sample"_corrected.fasta" $sample"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log

  
  source deactivate
}

