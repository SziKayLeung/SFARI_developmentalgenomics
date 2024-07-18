#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=vc362@exeter.ac.uk # email me at job completion
#SBATCH --mem=20G

export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/sequence/"

export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/"

export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/targeted/"


refAnno=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Reference/gencode.v38.annotation.gtf
refFile=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Reference/Hg38.fa

PROJECT="/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/split"

cd ${PROJECT}

#### Whole
#SQANTI3 QC

#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_chrX.gff" ${refAnno} ${refFile} -o Whole_cleaned_aligned_merged_collapsed_qced_chrX --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_chrY.gff" ${refAnno} ${refFile} -o Whole_cleaned_aligned_merged_collapsed_qced_chrY --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

#SQANTI3 filter
#Don't have faa file anymore as an input because we did --skiporf in the SqANTI3 QC step
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output Whole_cleaned_aligned_merged_collapsed_qced_filtered_chrX --skip_report --gtf "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_chrX_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_chrX_classification.txt"
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output Whole_cleaned_aligned_merged_collapsed_qced_filtered_chrY --skip_report --gtf "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_chrY_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_chrY_classification.txt"


#### Whole+Targeted

python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_chrX.gff" ${refAnno} ${refFile} -o WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrX --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_chrY.gff" ${refAnno} ${refFile} -o WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrY --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"


#SQANTI3 filter
#Don't have faa file anymore as an input because we did --skiporf in the SqANTI3 QC step
python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output WholeTargeted_cleaned_aligned_merged_collapsed_qced_filtered_chrX --skip_report --gtf "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrX_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrX_classification.txt"
python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output WholeTargeted_cleaned_aligned_merged_collapsed_qced_filtered_chrY --skip_report --gtf "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrY_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_chrY_classification.txt"

