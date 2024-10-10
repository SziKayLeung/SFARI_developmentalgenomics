#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion

# 08/10/2024: read stats of demux files

WHOLE=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/1_raw/PromethION/Project_10380/Upload_to_SRA
TARGETED=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined
OUTPUTDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_output

module load Miniconda2
source activate nanopore

# output_stats <input_dir> <output_file>
output_stats(){

	input_dir=$1
	output_file=$2

	# Clear the output file if it exists
	> "$output_file"

	# Loop through each FASTA file and run seqkit stats
	count=0
	for i in ${input_dir}/*fastq.gz*; do
		echo "Processing $i"
		
		if [ ${count} == 0 ]; then
			# Run seqkit stats and filter out the header line, appending results to the output file
			echo "first file"
			seqkit stats "$i" >> "$output_file"  # Skip the first line (header)
		else
			echo "next file"
			seqkit stats "$i" | awk 'NR > 1' >> "$output_file"  # Skip the first line (header)
		fi	
		(( count++ ))
		
	done
	
}

output_stats ${WHOLE} ${OUTPUTDIR}/WholeDemuxReadStats.txt
output_stats ${TARGETED} ${OUTPUTDIR}/TargetedDemuxReadStats.txt

#seqkit stats ${TARGETED}/*.fastq.gz -a > ${OUTPUTDIR}/TargetedDemuxReadStats.txt
