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
#SBATCH --array 235 # 1-383


export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/sequence/"

export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/"

export PYTHONPATH=$PYTHONPATH:"/lustre/home/vc362/lib/cDNA_Cupcake/targeted/"


refAnno=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Reference/gencode.v38.annotation.gtf
refFile=/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Reference/Hg38.fa

PROJECT="/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/split/"

cd ${PROJECT}

####  Targeted
#SQANTI3 QC
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Targeted_cleaned_aligned_merged_collapsed_chr${SLURM_ARRAY_TASK_ID}.gff" ${refAnno} ${refFile} -o Targeted_cleaned_aligned_merged_collapsed_qced_chr${SLURM_ARRAY_TASK_ID} --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"


#SQANTI3 filter
#Don't have faa file anymore as an input because we did --skiporf in the SqANTI3 QC step
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output Targeted_cleaned_aligned_merged_collapsed_qced_filtered_chr${SLURM_ARRAY_TASK_ID} --skip_report --gtf "${PROJECT}/Targeted_cleaned_aligned_merged_collapsed_qced_chr${SLURM_ARRAY_TASK_ID}_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/Targeted_cleaned_aligned_merged_collapsed_qced_chr${SLURM_ARRAY_TASK_ID}_classification.txt"

#### Whole
#SQANTI3 QC
#i=`sed -n "${SLURM_ARRAY_TASK_ID}p;d" /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/chunks.txt`
#grep -wF -f /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/${i} /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_${i::-12}.gff >> /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/Whole_cleaned_aligned_merged_collapsed_${i::-4}.gff

#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/Whole_cleaned_aligned_merged_collapsed_${i::-4}.gff" ${refAnno} ${refFile} -o Whole_cleaned_aligned_merged_collapsed_qced_${i::-4} --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

#SQANTI3 filter
#Don't have faa file anymore as an input because we did --skiporf in the SqANTI3 QC step
#python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output Whole_cleaned_aligned_merged_collapsed_qced_filtered_${i::-4} --skip_report --gtf "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_${i::-4}_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/Whole_cleaned_aligned_merged_collapsed_qced_${i::-4}_classification.txt"


#### Whole+Targeted
i=`sed -n "${SLURM_ARRAY_TASK_ID}p;d" /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/WholeTargeted_chunks.txt`
#grep -wF -f /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/${i} /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_${i:14:-12}.gff >> /gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/${i::-4}.gff

python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_qc.py "/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/split/${i::-4}.gff" ${refAnno} ${refFile} -o WholeTargeted_cleaned_aligned_merged_collapsed_qced_${i:14:-4} --report skip --genename --skipORF --CAGE_peak "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt"


#SQANTI3 filter
#Don't have faa file anymore as an input because we did --skiporf in the SqANTI3 QC step
python /lustre/home/vc362/lib/SQANTI3-5.1.2/sqanti3_filter.py rules --output WholeTargeted_cleaned_aligned_merged_collapsed_qced_filtered_${i:14:-4} --skip_report --gtf "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_${i:14:-4}_corrected.gtf" --json_filter "/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json" "${PROJECT}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_${i:14:-4}_classification.txt"

