#!/bin/bash

#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=12:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH --mem=5G # specify bytes memory to reserve
#SBATCH --array 1-19446:100


source /gpfs/ts0/shared/software/Miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate ficle

FICLE_ROOT=/lustre/home/vc362/lustre_project/ficle/FICLE/
export PATH=$PATH:${FICLE_ROOT}
export PATH=$PATH:${FICLE_ROOT}/reference

for s in $(seq $SLURM_ARRAY_TASK_ID `echo $SLURM_ARRAY_TASK_ID+99|bc`);
do

i=`sed -n "${s}p;d" /lustre/home/vc362/protein-coding-genes.txt` 
echo ${i}

cd /lustre/home/vc362/lustre_project/ficle/input_files

grep -w ${i} /lustre/recovered/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples.gtf.cds.gff > ${i}.gtf

/lustre/home/vc362/lustre_project/ficle/FICLE/reference/subset_reference_by_gene.py --r /lustre/home/vc362/resources/gencode.v44.annotation.gtf --g ${i} --o ./

gtfToGenePred ${i}.gtf ${i}.genePred
genePredToBed ${i}.genePred ${i}.bed12
sort -k1,1 -k2,2n ${i}.bed12 > ${i}_sorted.bed12

rm ${i}.bed12
rm ${i}.genePred

cd /lustre/home/vc362/lustre_project/ficle/

FICLE/ficle.py -r input_files/${i}_gencode.gtf -n ${i} -b input_files/${i}_sorted.bed12 -g input_files/${i}.gtf -c /lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_classification_2reads2samples_monomultirem.txt --cpat /lustre/recovered/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.ORF_prob.best.tsv -o results/
done
