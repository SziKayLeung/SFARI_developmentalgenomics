#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=80:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=cupcake_WholeTargeted_chr6_repeat2.o
#SBATCH --error=cupcake_WholeTargeted_chr6_repeat2.e

##-------------------------------------------------------------------------

echo Job started on:
date -u

module load Miniconda2/4.3.21

# paths
UTILS_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/0_utilities
MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom
TCLEAN_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/4_tclean
TCLEAN_WHOLE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/4_transcriptClean/Whole
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/6_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing

chrNum=chr6

##-------------------------------------------------------------------------

# merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y)
cd ${MERGED_CHROM_DIR}

#source activate isoseq3
#isoseq3 collapse ${MERGED_CHROM_DIR}/${chrNum}.bam ${chrNum}.gff --do-not-collapse-extra-5exons --min-aln-coverage=0.85 --min-aln-identity=0.95 --num-threads 16 --log-level=INFO --log-file=${chrNum}_collapse.log


for i in ${MERGED_CHROM_DIR}/chr6.read_stat.txt; do 
  chromosome=$(basename $i .read_stat.txt)
  echo $chromosome
  replaceONTprefix="${chromosome//chr/ONT}"
  echo $replaceONTprefix
    
  # replace gff with the correct gff name
  sed "s/PB/${replaceONTprefix}/g" $MERGED_CHROM_DIR/$chromosome.gff > $MERGED_CHROM_DIR/$chromosome.renamed.gff
   
  # replace read.stat txt with corrected isoform id
  sed "s/PB/${replaceONTprefix}/g" $MERGED_CHROM_DIR/$chromosome.read_stat.txt > $MERGED_CHROM_DIR/$chromosome.read_stat.renamed.txt

done


source activate nanopore
demux_cupcake_collapse.py ${MERGED_CHROM_DIR}/${chrNum}.read_stat.renamed.txt ${DEMUX_DIR}/WholeTargeted_sample_id.csv -o WholeTargeted_demux_${chrNum} --dataset ont -d ${DEMUX_DIR}