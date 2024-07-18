#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrchq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

# call orf
module load Miniconda2
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteogenomics_functions.sh
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteomics.config

mkdir -p $WKD_ROOT/5_calledOrfs; cd $WKD_ROOT/5_calledOrfs
source activate sqanti2_py3 


#python $LREAD/orf_calling/src/orf_calling.py \
#--orf_coord $NAME".ORF_prob.tsv" \
#--orf_fasta $NAME".ORF_seqs.fa" \
#--gencode $GENOME_GTF \
#--sample_gtf $ISO_GTF \
#--pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
#--classification $ISO_CLASSFILE \
#--sample_fasta $ISO_FASTA \
#--num_cores 2 \
#--output $NAME"_best_orf.tsv"


#mkdir -p $WKD_ROOT/6_refined_database; cd $WKD_ROOT/6_refined_database
#python $LREAD/refine_orf_database/src/refine_orf_database.py \
#  --name $NAME \
#  --orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv"  \
#  --pb_fasta $ISO_FASTA \
#  --coding_score_cutoff $coding_score_cutoff
 
  
#python $LREAD/visualization_track/src/make_pacbio_cds_gtf.py \
#--name $NAME \
#--sample_gtf $ISO_GTF \
#--refined_database $NAME"_orf_refined.tsv" \
#--called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
#--pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
#--include_transcript yes &> make_cds_gtf1.log

# modified make_pacbio_cds_gtf.py to remove cpm
#python $LREAD/visualization_track/src/make_pacbio_cds_gtf.py \
#--name $NAME"_no_transcript" \
#--sample_gtf $ISO_GTF \
#--refined_database $NAME"_orf_refined.tsv" \
#--called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
#--pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
#--include_transcript no &> make_cds_gtf1.log   
  
  
classify_protein
NAME=testing
ONAME=WholeTargeted
  
python $LREAD/rename_cds_to_exon/src/rename_cds_to_exon.py \
  --sample_gtf $WKD_ROOT/7_classified_protein/testing/$NAME"_with_cds.gtf" \
  --sample_name $NAME \
  --reference_gtf ./Trem2_Genome.gtf \
  --reference_name gencode \
  --num_cores 2 &> rename_cds_to_exon.log

python $LREAD/sqanti_protein/src/sqanti3_protein.py $NAME.transcript_exons_only.gtf \
  $NAME.cds_renamed_exon.gtf \
  $WKD_ROOT/5_calledOrfs/$ONAME"_best_orf.tsv" \
  gencode.transcript_exons_only.gtf \
  gencode.cds_renamed_exon.gtf \
  -d ./ \
  -p $NAME &> sqanti3_protein.log

# modified script with capture_output=TRUE as only works for python 3.7 and above
python $LREAD/5p_utr_status/src/1_get_gc_exon_and_5utr_info.py \
  --gencode_gtf ./Trem2_Genome.gtf \
  --odir ./

python $LREAD/5p_utr_status/src/2_classify_5utr_status.py \
  --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
  --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
  --sample_cds_gtf $WKD_ROOT/7_classified_protein/testing/$NAME"_with_cds.gtf" \
  --odir ./

python $LREAD/5p_utr_status/src/3_merge_5utr_info_to_pclass_table.py \
  --name $NAME \
  --utr_info pb_5utr_categories.tsv \
  --sqanti_protein_classification $NAME.sqanti_protein_classification.tsv \
  --odir ./

python $LREAD/protein_classification/src/protein_classification_add_meta.py \
  --protein_classification $NAME.sqanti_protein_classification_w_5utr_info.tsv \
  --best_orf $WKD_ROOT/5_calledOrfs/$ONAME"_best_orf.tsv" \
  --refined_meta $WKD_ROOT/6_refined_database/$ONAME"_orf_refined.tsv" \
  --ensg_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
  --name $NAME \
  --dest_dir ./

python $LREAD/protein_classification/src/protein_classification.py \
  --sqanti_protein $NAME.protein_classification_w_meta.tsv \
  --name $NAME"_unfiltered" \
  --dest_dir ./

python $LREAD/protein_gene_rename/src/protein_gene_rename.py \
  --sample_gtf $WKD_ROOT/7_classified_protein/testing/$NAME"_with_cds.gtf" \
  --sample_protein_fasta $WKD_ROOT/6_refined_database/$ONAME"_orf_refined.fasta" \
  --sample_refined_info $WKD_ROOT/6_refined_database/$ONAME"_orf_refined.tsv" \
  --pb_protein_genes $NAME"_genes.tsv" \
  --name $NAME

python $LREAD/protein_filter/src/protein_filter.py \
  --protein_classification $NAME"_unfiltered.protein_classification.tsv" \
  --gencode_gtf ./Trem2_Genome.gtf \
  --protein_fasta $NAME.protein_refined.fasta \
  --sample_cds_gtf $NAME"_with_cds_refined.gtf" \
  --min_junctions_after_stop_codon $min_junctions_after_stop_codon \
  --name $NAME
