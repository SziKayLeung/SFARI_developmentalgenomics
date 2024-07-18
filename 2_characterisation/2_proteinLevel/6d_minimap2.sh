#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=80:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=6d_minimap2.o
#SBATCH --error=6d_minimap2.e

# 01/12/2023: Minimap alignment of best orf to genome fasta for isoform visulisation

#-----------------------------------------------------------------------#

# load and source packages
module load Miniconda2

source activate sqanti2_py3
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteogenomics_functions.sh
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics/1_characterisation/2_proteinLevel/0_proteomics.config

LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PYTHONPATH=$PYTHONPATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/cDNA_Cupcake/sequence
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/cDNA_Cupcake/sequence
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/cDNA_Cupcake

# extract best orf
source activate sqanti2_py3  
cd $WKD_ROOT/5_calledOrfs
#awk '{print $2}' $WKD_ROOT/5_calledOrfs/$NAME".ORF_prob.best.tsv" > $WKD_ROOT/5_calledOrfs/$NAME".ORF_prob.best_ID.tsv"
#seqtk subseq $WKD_ROOT/5_calledOrfs/$NAME".ORF_seqs.fa" $WKD_ROOT/5_calledOrfs/$NAME".ORF_prob.best_ID.tsv" > $WKD_ROOT/5_calledOrfs/$NAME"_bestORF".fa
##extract_fasta_bestorf.py --fa $WKD_ROOT/5_calledOrfs/$NAME".ORF_seqs.fa" --orf $WKD_ROOT/5_calledOrfs/$NAME".ORF_prob.best.tsv" --o_name $NAME"_bestORF" --o_dir $WKD_ROOT/5_calledOrfs 


# write.table(proteinInput$t.class.files[proteinInput$t.class.files$base_acc == "ONT20_1826_53630","isoform"],"/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/6a_longReadProteogenomics/5_calledOrfs/ONT20_1826_53630_collapsed_ID.txt", quote=F, row.names=F)
#grep -f ONT20_1826_53630_collapsed_ID.txt Whole.ORF_prob.best.tsv | awk '{print $2}' > ONT20_1826_53630_ORF_prob.best_ID.tsv
#seqtk subseq $WKD_ROOT/5_calledOrfs/$NAME"_bestORF".fa ONT20_1826_53630_ORF_prob.best_ID.tsv > $WKD_ROOT/5_calledOrfs/ONT20_1826_53630_collapsed_ID.fa
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA ONT20_1826_53630_collapsed_ID.fa > ONT20_1826_53630.sam 2> ONT20_1826_53630.map.log
samtools sort -O SAM ONT20_1826_53630.sam > ONT20_1826_53630.sorted.sam
sam_to_gff3.py ONT20_1826_53630.sam -s ${SPECIES}
gffread ONT20_1826_53630.gff3 -T -o ONT20_1826_53630.gtf

echo "All reads"
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $NAME"_bestORF.fa" > $NAME.sam 2> $NAME.map.log
samtools sort -O SAM $NAME.sam > $NAME.sorted.sam
sam_to_gff3.py $NAME.sam -s ${SPECIES}
gffread $NAME.gff3 -T -o $NAME.gtf