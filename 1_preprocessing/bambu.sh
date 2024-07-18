#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --output=runBambu.o
#SBATCH --error=runBambu.e


# 18/07/2024: Bambu, sqanti3 from pbmm2 aligned files onwards 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen

# input variables 
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.gtf
CAGE="/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed"
POLYA="/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
alignedDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/5_isoseq/pbmm2_align
outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/16_bambu
CUPCAKE_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake
SQANTI3_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/software/SQANTI3
SQANTIjson="/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" 

##-------------------------------------------------------------------------

module load R/4.2.2-foss-2022b
Rscript ${LOGEN_ROOT}/assist_ont_processing/run_Bambu.R -i ${alignedDir} -f ${GENOME_FASTA} -g ${GENOME_GTF} -o ${outputDir} --index "aligned.bam"

source activate sqanti2_py3
cd ${outputDir}
mkdir -p sqanti; cd sqanti
sample=WholeTargeted
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_DIR/sequence

# In case this is still an issue for anyone, this code removes anything that is not + or - in the strand column:
awk '$7 != "." || NR < 3' ${outputDir}/extended_annotations.gtf > ${outputDir}/extended_corrected_annotations.gtf
python $SQANTI3_DIR/sqanti3_qc.py ${outputDir}/extended_corrected_annotations.gtf ${GENOME_GTF} ${GENOME_FASTA} -o $sample --report skip --genename --skipORF --CAGE_peak ${CAGE} --polyA_motif_list ${POLYA}
python $SQANTI3_DIR/sqanti3_filter.py rules $sample"_classification.txt" --isoforms $sample"_corrected.fasta" --gtf $sample"_corrected.gtf" --json_filter ${SQANTIjson}