#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --mem=20G
#SBATCH --array=1-343%20
#SBATCH --output=log_sqanti/sqanti-%A_%a.o
#SBATCH --error=log_sqanti/sqanti-%A_%a.e

# 04/09/2024: running sqanti from split files

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
source activate sqanti2_py3

SPLIT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/split
SQANTI_OUTPUT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti
MERGED_CHROM_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/WholeTargeted/mergedChrom

## Software 
export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export CUPCAKE=${SOFTDIR}/cDNA_Cupcake
export SEQUENCE=$CUPCAKE/sequence
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export SQANTI3_DATA_DIR=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data
export CAGE_PEAK=${SQANTI3_DATA_DIR}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=${SQANTI3_DATA_DIR}/polyA_motifs/mouse_and_human.polyA_motif.txt
export refAnno=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/gencode.v38.annotation.gtf
export refFile=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/Hg38.fa


##-------------------------------------------------------------------------

cd ${SQANTI_OUTPUT}

# print the SLURM_ARRAY_TASK_ID nth line and delete the rest
i=`sed -n "${SLURM_ARRAY_TASK_ID}p;d" ${SPLIT}/WholeTargeted_chunks.txt`
prefix=$(basename $i .txt)

## extract variables from $i for downstream labelling 
## e.g. i=WholeTargeted_chr10_chunk00.txt
## chr="chr10_chunk00.txt" (extract after WholeTargeted)
chr=${i#WholeTargeted_}

## chr_chunk=chunk10_chunk00 (remove txt before chr10_chunk00.txt)
chr_chunk=${chr%.txt}

## remove everything after (and including) the first underscore _ in the string of chr10_chunk00.txt
## chr=chr10
chr=${chr%%_*}

echo "Processing $i; output chromosome = $chr; output gff: $prefix.gff; output sqanti: $chr_chunk"
grep -wF -f ${SPLIT}/${i} ${MERGED_CHROM_DIR}/${chr}.gff >> ${SPLIT}/${prefix}.gff

python ${SQANTI3_DIR}/sqanti3_qc.py ${SPLIT}/${prefix}.gff ${refAnno} ${refFile} -o WholeTargeted_collapsed${chr_chunk} \
  --report skip --genename --skipORF \
  --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA}