#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=144:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=am1248@exeter.ac.uk # email me at job completion
#SBATCH --output=pancreasQC.o
# #SBATCH --job-name=squanti3-%A_%a_%u
#SBATCH --mem=200G

module load Miniconda2
source activate sqanti2_py3

export PYTHONPATH=$PYTHONPATH:/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Resources/SQANTI3/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Resources/SQANTI3/cDNA_Cupcake/
##export PYTHONPATH=$PYTHONPATH:"/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-4.2/cDNA_Cupcake/targeted/"
export SOFTDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software
export SQANTI3_DIR=$SOFTDIR/SQANTI3
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
  
refAnno=/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Resources/gencode.v41.annotation.gtf
refFile=/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Resources/Hg38.fa
gff=/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/merged_collapse.gff
counts=/gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/demux_fl_count.csv
cage=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
polyA=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt

cd /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/
python $SQANTI3_DIR/sqanti3_qc.py merged_collapse.gff ${refAnno} ${refFile} -o final --genename --fl_count ${counts} --CAGE_peak ${cage} --polyA_motif_list ${polyA}

python /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/SQANTI/SQANTI3-5.1.1/sqanti3_filter.py rules \
--output SQANTI3_whole \
--faa /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/final_corrected.faa \
--gtf /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/final_corrected.gtf \
--json_filter /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Resources/SQANTI3-5.1/utilities/filter/filter_default.json /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/final_classification.txt

source activate nanopore
subset_quantify_filter_tgenes.R  \
--classfile /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/SQANTI3_whole_RulesFilter_result_classification.txt \
--expression /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/demux_fl_count.csv \
--target_genes /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Target2.txt

subset_quantify_filter_tgenes.R  \
--classfile /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/SQANTI3_whole_RulesFilter_result_classification.txt \
--expression /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Ailsa/AM/20220623_1515_2E_PAI82330_18b0e941/Porechop/pc_test/combined/demux_fl_count.csv \
--target_genes /gpfs/mrc0/projects/Research_Project-MRC190311/Ailsa/Targeted/Target2.txt \
--filter --nsample=1 --nreads=2