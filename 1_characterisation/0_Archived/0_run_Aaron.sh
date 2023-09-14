#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=0_run_Aaron.o
#SBATCH --error=0_run_Aaron.e

# 04/07/2023: Run Aaron's file for testing pipeline

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
SFARI_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/SFARI_developmentalgenomics
source $SFARI_ROOT/1_characterisation/01_source_functions.sh
source $SFARI_ROOT/1_characterisation/SFARI_characterisation.config
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing


##-------------------------------------------------------------------------

export dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Collabs/aaronOptimisation
export sample=pyc
rawdata=/gpfs/mrc0/projects/Research_Project-MRC148213/Aaron/ONT/pyc.zip
#cp ${rawdata} ${dir}
#unzip ${dir}/${rawdata}

export fastq=${dir}/pychopper/full_length_output.fq

mkdir -p ${dir}/1_align ${dir}/2_collapse ${dir}/3_sqanti3

##-------------------------------------------------------------------------
# Run pipeline

source activate isoseq3

# align
#cd ${dir}/1_align
#pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${fastq} ${sample}_mapped.bam --unmapped --log-level TRACE --log-file ${sample}_mapped.log

# collapse
#cd ${dir}/2_collapse
#isoseq3 collapse ${dir}/1_align/${sample}_mapped.bam ${sample}.gff --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons --log-level TRACE --log-file ${sample}_collapsed.log

source activate sqanti2_py3

# sqanti3
cd ${dir}/3_sqanti3
#python $SQANTI3_DIR/sqanti3_qc.py -t 30 ${dir}/2_collapse/${sample}.gff $GENOME_GTF $GENOME_FASTA --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip &> ${sample}.sqanti.qc.log

# sqanti3 filter
python $SQANTI3_DIR/sqanti3_filter.py rules ${sample}"_classification.txt" \
--faa=${sample}"_corrected.fasta" \
--gtf=${sample}"_corrected.gtf" \
-j=${reducedJson}  &> ${sample}.sqanti.filter.log

sqUtilDir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/SQANTI3/utilities
Rscript ${sqUtilDir}/report_qc/SQANTI3_report.R ${dir}/3_sqanti3/${sample}_RulesFilter_result_classification.txt ${dir}/3_sqanti3/${sample}_junctions.txt ${sqUtilDir} False pdf
