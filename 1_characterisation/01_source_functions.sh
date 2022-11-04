
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
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report pdf &> $sample.sqanti.qc.log
  
  echo "Processing Sample $sample for SQANTI2 filter"
  python $SQANTI3_DIR/sqanti3_RulesFilter.py $sample"_classification.txt" $sample"_corrected.fasta" $sample"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log

  
  source deactivate
}


# subset_case_control <class_file> <meta_file> <abundance_file> <output_dir> <cpat_file> <cpat_noORF> <control_variable> <case_variable>
subset_case_control(){
  
  # variables 
  class_file=$1
  meta_file=$2
  abundance_file=$3
  output_dir=$4
  cpat_file=$5
  cpat_noORF=$6
  control=$7
  case=$8
  
  source activate nanopore
  cond=(Case Control Case_Control)
  
  ### Subset by condition and by sample
  # subset cases and control by counts
  Rscript $SQTABCOUNTS -f ${class_file} -m ${meta_file} -c ${abundance_file} -o ${output_dir}
  class_count_file=$(ls ${output_dir}/*classification.filtered_lite_classification_counts.txt)
  echo "********* Subsetting condition by counts"
  echo "Using: $class_count_file"
  Rscript $SQCOUNT -f ${class_count_file} -m ${meta_file} -o ${output_dir} -d partial
  
  ### generate classification.txt, gtf file, and corresponding SQANTI report
  sq_dir=$(dirname ${class_file})
  sqname=$(basename ${class_file})
  sqname_prefix=${sqname//"_classification.txt"/}
  
  # for each condition (case, control, case_control)
  for c in ${cond[@]}; do 
  
    echo "Processing $c"
    
    Rscript $SQSUBSET -i ${output_dir}/$c"_ID.txt" -d $sq_dir -s $sqname_prefix -n $c -o ${output_dir} -f 1
    Rscript $SQ_Report ${output_dir}/$c"_classification.txt" ${output_dir}/$c"_junctions.txt"
    
    # convert gtf to bed file 
    convert_gtf_bed12 ${output_dir}/$c.gtf 
    
    # colour the bed file by abundance   
    source activate nanopore  
    if [ $c == "Case_Control" ] ; then
      python $ISOCOL --bed ${output_dir}/$c"_sorted.bed12" --cpat ${cpat_file} --noORF ${cpat_noORF} --a ${output_dir}/$c"_CaseAbundance.csv" --o ${output_dir}/$c"_caseabundance_sorted_coloured.bed12"
      python $ISOCOL --bed ${output_dir}/$c"_sorted.bed12" --cpat ${cpat_file} --noORF ${cpat_noORF} --a ${output_dir}/$c"_ControlAbundance.csv" --o ${output_dir}/$c"_controlabundance_sorted_coloured.bed12"
      else
        python $ISOCOL --bed ${output_dir}/$c"_sorted.bed12" --cpat ${cpat_file} --noORF ${cpat_noORF} --a ${output_dir}/$c"_Only_Abundance.csv" --o ${output_dir}/$c"_sorted_coloured.bed12"
    fi
  
  done
  
  echo "Control = " ${control} > ${output_dir}/README.md
  echo "Case =" ${case} >> ${output_dir}/README.md
  
}



