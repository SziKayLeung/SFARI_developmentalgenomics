#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=reOrientateTargeted22.o
#SBATCH --error=reOrientateTargeted22.e


# 09/06/2024: reorientate sample 22 (duplicated with targeted06)

##-------------------------------------------------------------------------

# source function
module load Miniconda2/4.3.21
source activate nanopore

WKD_ROOT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome
RAW_FASTQ=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Targeted_transcriptome/UploadtoSRA/combined
SUBSETPOLYTAILS=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/assist_ont_processing/subset_polyA_polyT.py

# 4) post_porechop_run_cutadapt <input_fastq> <output_dir>
post_porechop_run_cutadapt(){
  
  input_dir=$(dirname $1)
  name=$(basename $1 .fastq.gz)
  echo "Output to $name"
  
  # requires fasta files for downstream
  echo "Converting $1 to fasta"
  gunzip -c $1 | seqtk seq -a - > ${input_dir}/${name}.fasta
  
  # subset fasta file to polyA and polyT fasta (i.e. reads ending with PolyA and starting with polyT)
  # reads that end with AAAAAAAAAA = plus reads 
  # reads that start with TTTTTTTTTT = minus reads (need to be reverse complemented)
  echo "Subsetting fasta to polyA and polyT sequences"
  python ${SUBSETPOLYTAILS} --fa ${input_dir}/${name}.fasta --o_name ${name} --o_dir $2
  
  # working in output directory
  cd $2
  
  # reverse complement minus reads (reads ending with polyT)
  seqtk seq -r ${name}_PolyT.fasta > ${name}_PolyT_rev.fasta
  
  # use cutadapt package to trim polyA
  echo "Remove polyA sequences using cutadapt"
  cutadapt -a "A{60}" ${name}_PolyA.fasta -o ${name}_PolyA_cutadapted.fasta &> ${name}_polyA_cutadapt.log
  cutadapt -a "A{60}" ${name}_PolyT_rev.fasta -o ${name}_PolyT_rev_cuptadapted.fasta &> ${name}_polyT_cutadapt.log
  
  # concatenated reverse minus polyT and polyA reads
  cat ${name}_PolyA_cutadapted.fasta ${name}_PolyT_rev_cuptadapted.fasta > ${name}_combined.fasta
  
}


##-------------------------------------------------------------------------

SamplePath=${RAW_FASTQ}/Targeted22.fastq.gz
Sample=$(basename ${SamplePath} .fastq.gz)
echo "Processing ${Sample}"

# delinate polyA and polyT sequences, reverse complement polyT sequences, remove polyA from all sequences
post_porechop_run_cutadapt ${SamplePath} ${WKD_ROOT}/2_cutadapt_merge

grep "^>" ${WKD_ROOT}/2_cutadapt_merge/${Sample}_combined.fasta | awk '{print $1}' | sed 's/>//g' > ${WKD_ROOT}/2_cutadapt_merge/${Sample}_combined_Ids.txt


