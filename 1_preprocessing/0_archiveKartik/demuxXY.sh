#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=120:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=v.chundru@exeter.ac.uk # email address
#SBATCH --array 2

#if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
#python /lustre/home/vc362/MRC148213/vc362/LOGen/assist_ont_processing/demux_cupcake_collapse.py /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_chrX.read_stat.renamed.txt /lustre/home/vc362/MRC148213/Rosie/WholeTargeted//BAMford/WholeTargeted_sample_id_fixed.csv -o WholeTargeted_demux_chrX
#else 
#python /lustre/home/vc362/MRC148213/vc362/LOGen/assist_ont_processing/demux_cupcake_collapse.py /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/WholeTargeted_cleaned_aligned_merged_collapsed_chrY.read_stat.renamed.txt /lustre/home/vc362/MRC148213/Rosie/WholeTargeted//BAMford/WholeTargeted_sample_id_fixed.csv -o WholeTargeted_demux_chrY
#fi

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
python /lustre/home/vc362/MRC148213/vc362/LOGen/assist_ont_processing/demux_cupcake_collapse.py /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_chrX.read_stat.renamed.txt /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/BAMford/Whole_sample_id_fixed.csv -o Whole_demux_chrX
else 
python /lustre/home/vc362/MRC148213/vc362/LOGen/assist_ont_processing/demux_cupcake_collapse.py /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/cleaned_merged_collapsed/Whole_cleaned_aligned_merged_collapsed_chrY.read_stat.renamed.txt /lustre/home/vc362/MRC148213/Rosie/WholeTargeted/BAMford/Whole_sample_id_fixed.csv -o Whole_demux_chrY
fi
