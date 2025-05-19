#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=60:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=150G # specify bytes memory to reserve
#SBATCH --array 23 # 24 chromsomes, 22 autosomal, X and Y
#SBATCH --output=../log/log_Oct2024/6_demux-%A_%a.o
#SBATCH --error=../log/log_Oct2024/6_demux-%A_%a.e
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

##-------------------------------------------------------------------------

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
DEMUX_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/8_demux

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing
mkdir -p ${DEMUX_DIR}/split

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21
source activate nanopore

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
chrNum=${chromosomes[${SLURM_ARRAY_TASK_ID}]}  

if compgen -G "${ISOSEQ_COLLAPSE_DIR}/${chrNum}.read_stat.renamed.txt" > /dev/null; then

    # Extract the header
    header=$(head -n 1 ${ISOSEQ_COLLAPSE_DIR}/${chrNum}.read_stat.renamed.txt)
    
    # Split the file into chunks, skipping the header for the split but keeping it
    tail -n +2 ${ISOSEQ_COLLAPSE_DIR}/${chrNum}.read_stat.renamed.txt | \
    split -l100000 - ${DEMUX_DIR}/split/${chrNum}.read_stat.renamed_chunk --additional-suffix=.txt -d

    # Add the header to each of the split files
    for file in ${DEMUX_DIR}/split/${chrNum}.read_stat.renamed_chunk*.txt; do
        (echo "$header" && cat "$file") > temp && mv temp "$file"
    done

    # List the split files
    ls ${DEMUX_DIR}/split/*${chrNum}.read_stat.renamed_chunk* > ${DEMUX_DIR}/split/${chrNum}_Allchunks_list.txt

    if compgen -G "$DEMUX_DIR/*${chrNum}_fl*.csv" > /dev/null; then
        echo WholeTargeted_demux_${chrNum}_fl_count.csv present
    else
        echo Run demux for $chrNum

        while read -r file; do
            prefix=$(basename $file .txt)
            echo $prefix

            demux_cupcake_collapse.py ${file} ${DEMUX_DIR}/WholeTargeted_sample_id.csv -o ${prefix} --dataset ont -d ${DEMUX_DIR}/split
       done < ${DEMUX_DIR}/split/${chrNum}_Allchunks_list.txt
   fi
else
    echo "$chrNum still collapsing"
fi
