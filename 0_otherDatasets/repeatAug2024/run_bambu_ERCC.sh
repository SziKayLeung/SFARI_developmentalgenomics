#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per nodee
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --output=2_runBambu.o
#SBATCH --error=2_runBambu.e

# Aim: run bambu in ERCC dataset

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen

# input variables 
alignedDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/repeatAug2024/ERCC/5_cupcake/5_align
outputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/repeatAug2024/ERCC/5_bambu
config=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LRPipeline/ont_cDNA/config/config_ERCC.yaml

##-------------------------------------------------------------------------

module load R/4.2.2-foss-2022b

source ${config}
Rscript ${LOGEN_ROOT}/assist_ont_processing/run_Bambu.R -i ${alignedDir} -g ${GENOME_GTF} --index="filtered.bam" -f ${GENOME_FASTA} -o ${outputDir} -q="FALSE"

export GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/ERCC/ERCC92.fa
export GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/ERCC/ERCC92.gtf
annotations <- prepareAnnotations("/lustre/projects/Research_Project-MRC148213/lsl693/references/ERCC/ERCC92.gtf")
test.bam <- list.files(path = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/repeatAug2024/ERCC/5_cupcake/5_align", 
pattern = paste0("filtered.bam", "$"), full = T)

se <- bambu(reads = test.bam, annotations = annotations, genome = "/lustre/projects/Research_Project-MRC148213/lsl693/references/ERCC/ERCC92.fa", quant = FALSE)
