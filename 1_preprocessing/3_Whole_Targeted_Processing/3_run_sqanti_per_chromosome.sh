#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=1# specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array 0-23 # 24 chromsomes, 22 autosomal, X and Y
#SBATCH --output=../log/log_Oct2024/3_run_sqanti_per_chromosome-%A_%a.o
#SBATCH --error=../log/log_Oct2024/3_run_sqanti_per_chromosome-%A_%a.e

# 17/10/2024: modify script to run sqanti through chromosome 

##-------------------------------------------------------------------------

module load Miniconda2/4.3.21

ISOSEQ_COLLAPSE_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/6_isoseqCollapse
SPLIT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/split
SQANTI_OUTPUT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax

## Software 
export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export CUPCAKE=${SOFTDIR}/cDNA_Cupcake
export SEQUENCE=$CUPCAKE/sequence
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export SQANTI3_DATA_DIR=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/data
export CAGE_PEAK=${SQANTI3_DATA_DIR}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=${SQANTI3_DATA_DIR}/polyA_motifs/mouse_and_human.polyA_motif.txt
export refAnno=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/gencode.v38.annotation.gtf
export refFile=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/Reference/Hg38.fa
export SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
chrNum=${chromosomes[${SLURM_ARRAY_TASK_ID}]}  

# modify_files <read_stat.txt> <output_split_directory>
modify_files(){
  chromosome=$(basename $1 .read_stat.txt)
  replaceONTprefix="${chromosome//chr/ONT}"
  echo $replaceONTprefix
  
  cd $2
  
  # extract the isoform list from each read stat file
  awk -F'\t' '{print $2}' $1 | tail -n +2 > WholeTargeted_${chromosome}.id.txt
  
  # replace PB with ONT and the chromosome 
  sed -i "s/PB/${replaceONTprefix}/g" WholeTargeted_${chromosome}.id.txt
  
  # replace gff with the correct gff name
  sed "s/PB/${replaceONTprefix}/g" $ISOSEQ_COLLAPSE_DIR/${chromosome}.gff > $ISOSEQ_COLLAPSE_DIR/${chromosome}.renamed.gff
  
  # replace read.stat txt with corrected isoform id
  sed "s/PB/${replaceONTprefix}/g" $ISOSEQ_COLLAPSE_DIR/${chromosome}.read_stat.txt > $ISOSEQ_COLLAPSE_DIR/${chromosome}.read_stat.renamed.txt
  
  # split chunk with 1000000 lines
  split WholeTargeted_${chromosome}.id.txt WholeTargeted_${chromosome}_chunk -l1000000 --additional-suffix=.txt -d
  
}

# run_sqanti_per_chunk <chunk_file> <input_collapse_dir> <input_split_dir> <output_sqanti_dir>
run_sqanti_per_chunk(){

	file=$(basename $1)
	prefix=$(basename $file .txt)

	## extract variables from $i for downstream labelling 
	## e.g. i=WholeTargeted_chr10_chunk00.txt
	## chr="chr10_chunk00.txt" (extract after WholeTargeted)
	chr=${file#WholeTargeted_}

	## chr_chunk=chunk10_chunk00 (remove txt before chr10_chunk00.txt)
	chr_chunk=${chr%.txt}

	## remove everything after (and including) the first underscore _ in the string of chr10_chunk00.txt
	## chr=chr10
	chr=${chr%%_*}

	echo "Processing $file; output chromosome = $chr; output gff: $prefix.gff; output sqanti: $chr_chunk"
	
	cd $4 
	## extract the relevant gff to run sqanti
	## renamed.gff from renaming the PB to ONT<chr> in the gff after cupcake collapse
	grep -wF -f $1 $2/${chr}.renamed.gff >> $3/${prefix}.gff

	## run sqanti
	echo "Run Sqanti"
	source activate sqanti2_py3
	python ${SQANTI3_DIR}/sqanti3_qc.py $3/${prefix}.gff ${refAnno} ${refFile} -o WholeTargeted_collapsed${chr_chunk} \
	  --report skip --genename \
	  --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA}
	  
	## run sqanti filter
	echo "Run Sqanti filter"
	python ${SQANTI3_DIR}/sqanti3_filter.py rules \
	  --output WholeTargeted_collapsed_filtered${chr_chunk} \
	  --skip_report --gtf WholeTargeted_collapsed${chr_chunk}_corrected.gtf \
	  --json_filter ${SQANTI_JSON} WholeTargeted_collapsed${chr_chunk}_classification.txt
}


if [ -f $ISOSEQ_COLLAPSE_DIR/${chrNum}.read_stat.txt ]; then
	if [ -f $ISOSEQ_COLLAPSE_DIR/${chrNum}.read_stat.renamed.txt ]; then
		echo "${chrNum} files already modified"
	else 
		echo "${chrNum} to modify files"
		modify_files $ISOSEQ_COLLAPSE_DIR/${chrNum}.read_stat.txt ${SPLIT}
		ls ${SPLIT}/*WholeTargeted_${chrNum}_chunk* > ${SPLIT}/WholeTargeted_${chrNum}_Allchunks_list.txt
	fi

else
	
	echo "${chrNum} still collapsing"
	exit 2

fi

while read -r file; do

        filename=$(basename $file .txt)
        chr_chunk=${filename#WholeTargeted_}

        if [ -f ${SQANTI_OUTPUT}/*${chr_chunk}_RulesFilter_result_classification.txt* ]; then
                echo "No need to run sqanti on $chr_chunk as already processed"
        else
                echo "Run sqanti on ${chr_chunk}"
                echo "To process $file"
                run_sqanti_per_chunk ${file} ${ISOSEQ_COLLAPSE_DIR} ${SPLIT} ${SQANTI_OUTPUT}
        fi

done < ${SPLIT}/WholeTargeted_${chrNum}_Allchunks_list.txt



