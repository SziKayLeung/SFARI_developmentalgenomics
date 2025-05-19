# prepare reference tables from reference files in config file
# generate files for downstream usage
prepare_reference_tables(){
  mkdir -p $WKD_ROOT/3_reference_tables; cd $WKD_ROOT/3_reference_tables
  
  echo "#*************************************** Preparing reference tables"
  python $LREAD/generate_reference_tables/src/prepare_reference_tables.py \
  --gtf $GENOME_GTF \
  --fa $GENOME_TRANSCRIPT_FASTA \
  --ensg_gene ensg_gene.tsv   \
  --enst_isoname enst_isoname.tsv \
  --gene_ensp gene_ensp.tsv \
  --gene_isoname gene_isoname.tsv \
  --isoname_lens isoname_lens.tsv \
  --gene_lens gene_lens.tsv \
  --protein_coding_genes protein_coding_genes.txt
  
  echo "Make gencode databases"
  python $LREAD/make_gencode_database/src/make_gencode_database.py \
  --gencode_fasta $GENOME_TRANSLATION_FASTA \
  --protein_coding_genes protein_coding_genes.txt \
  --output_fasta gencode_protein.fasta \
  --output_cluster gencode_isoname_clusters.tsv
}

summarise_longread_data(){
  mkdir -p $WKD_ROOT/4_longread; cd $WKD_ROOT/4_longread
  
  source activate sqanti2_py3
  echo "#*************************************** Generate fasta file with six nucleotide frame translatation"
  python $LREAD/six_frame_translation/src/six_frame_translation.py \
  --iso_annot $ISO_CLASSFILE \
  --ensg_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
  --sample_fasta $ISO_FASTA \
  --output_fasta $NAME.6frame.fasta
  
  echo "#*************************************** Summarise long-read transcriptome"
  python $LOGEN/transcriptome_summary.py \
  --sq_out $ISO_CLASSFILE \
  --ensg_to_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
  --enst_to_isoname $WKD_ROOT/3_reference_tables/enst_isoname.tsv \
  --len_stats $WKD_ROOT/3_reference_tables/gene_lens.tsv

}

call_orf(){
  echo "#*************************************** Calling open reading frames from long-read transcriptome data"
  mkdir -p $WKD_ROOT/5_calledOrfs; cd $WKD_ROOT/5_calledOrfs
  
  source activate sqanti2_py3 
  cpat.py \
  -x $HEXAMER \
  -d $LOGITMODEL \
  -g $ISO_FASTA \
  --min-orf=50 \
  --top-orf=50 \
  -o $NAME \
  1> $NAME"_cpat.output" \
  2> $NAME"_cpat.error"
  
  python $LREAD/orf_calling/src/orf_calling.py \
  --orf_coord $NAME".ORF_prob.tsv" \
  --orf_fasta $NAME".ORF_seqs.fa" \
  --gencode $GENOME_GTF \
  --sample_gtf $ISO_GTF \
  --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
  --classification $ISO_CLASSFILE \
  --sample_fasta $ISO_FASTA \
  --num_cores 2 \
  --output $NAME"_best_orf.tsv"
  
  source deactivate 

}

refine_calledorf(){
  echo "#*************************************** Filter on called open reading frames"
  mkdir -p $WKD_ROOT/6_refined_database; cd $WKD_ROOT/6_refined_database
  
  source activate sqanti2_py3
  python $LOGEN/refine_orf_database.py \
  --name $NAME \
  --orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv"  \
  --pb_fasta $ISO_FASTA \
  --coding_score_cutoff $coding_score_cutoff &> refine_org.log
  
  python $LOGEN/make_pacbio_cds_gtf.py \
  --name $NAME \
  --sample_gtf $ISO_GTF \
  --refined_database $NAME"_orf_refined.tsv" \
  --called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
  --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
  --include_transcript yes &> make_cds_gtf1.log
  
  # modified make_pacbio_cds_gtf.py to remove cpm
  python $LOGEN/make_pacbio_cds_gtf.py \
  --name $NAME"_no_transcript" \
  --sample_gtf $ISO_GTF \
  --refined_database $NAME"_orf_refined.tsv" \
  --called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
  --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
  --include_transcript no &> make_cds_gtf1.log   
}


classify_protein(){
  mkdir -p $WKD_ROOT/7_classified_protein; cd $WKD_ROOT/7_classified_protein
  
  source activate sqanti2_py3
  
  python $LREAD/rename_cds_to_exon/src/rename_cds_to_exon.py \
  --sample_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
  --sample_name $NAME \
  --reference_gtf $GENOME_GTF \
  --reference_name gencode \
  --num_cores 2 &> rename_cds_to_exon.log
  
  python $LOGEN/sqanti3_protein.py $NAME.transcript_exons_only.gtf \
  $NAME.cds_renamed_exon.gtf \
  $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
  gencode.transcript_exons_only.gtf \
  gencode.cds_renamed_exon.gtf \
  -d ./ \
  -p $NAME &> sqanti3_protein.log
  
  # modified script with capture_output=TRUE as only works for python 3.7 and above
  python $LREAD/5p_utr_status/src/1_get_gc_exon_and_5utr_info.py \
  --gencode_gtf $GENOME_GTF \
  --odir ./
    
    python $LREAD/5p_utr_status/src/2_classify_5utr_status.py \
  --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
  --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
  --sample_cds_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
  --odir ./
    
    python $LREAD/5p_utr_status/src/3_merge_5utr_info_to_pclass_table.py \
  --name $NAME \
  --utr_info pb_5utr_categories.tsv \
  --sqanti_protein_classification $NAME.sqanti_protein_classification.tsv \
  --odir ./
    
    python $LREAD/protein_classification/src/protein_classification_add_meta.py \
  --protein_classification $NAME.sqanti_protein_classification_w_5utr_info.tsv \
  --best_orf $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
  --refined_meta $WKD_ROOT/6_refined_database/$NAME"_orf_refined.tsv" \
  --ensg_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
  --name $NAME \
  --dest_dir ./
    
    python $LREAD/protein_classification/src/protein_classification.py \
  --sqanti_protein $NAME.protein_classification_w_meta.tsv \
  --name $NAME"_unfiltered" \
  --dest_dir ./
    
    python $LREAD/protein_gene_rename/src/protein_gene_rename.py \
  --sample_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
  --sample_protein_fasta $WKD_ROOT/6_refined_database/$NAME"_orf_refined.fasta" \
  --sample_refined_info $WKD_ROOT/6_refined_database/$NAME"_orf_refined.tsv" \
  --pb_protein_genes $NAME"_genes.tsv" \
  --name $NAME
  
  python $LREAD/protein_filter/src/protein_filter.py \
  --protein_classification $NAME"_unfiltered.protein_classification.tsv" \
  --gencode_gtf $GENOME_GTF \
  --protein_fasta $NAME.protein_refined.fasta \
  --sample_cds_gtf $NAME"_with_cds_refined.gtf" \
  --min_junctions_after_stop_codon $min_junctions_after_stop_codon \
  --name $NAME
}

run_hybrid_annotation(){
  mkdir -p $WKD_ROOT/8_hybrid_annotation; cd $WKD_ROOT/8_hybrid_annotation
  python $LOGEN/make_hybrid_database.py \
  --protein_classification $WKD_ROOT/7_classified_protein/$NAME".classification_filtered.tsv" \
  --gene_lens $WKD_ROOT/3_reference_tables/gene_lens.tsv \
  --pb_fasta $WKD_ROOT/7_classified_protein/$NAME".filtered_protein.fasta" \
  --gc_fasta $WKD_ROOT/3_reference_tables/gencode_protein.fasta \
  --refined_info $WKD_ROOT/7_classified_protein/$NAME"_orf_refined_gene_update.tsv" \
  --pb_cds_gtf $WKD_ROOT/7_classified_protein/$NAME"_with_cds_filtered.gtf" \
  --name $NAME \
  --lower_kb $lower_kb \
  --upper_kb $upper_kb \
  --lower_cpm $lower_cpm
}
