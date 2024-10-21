#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=log/merged_post_sqanti.o
#SBATCH --error=log/merged_post_sqanti.e


# 09/09/2024: post sqanti per chromosome

##-------------------------------------------------------------------------

# paths and variables
SQANTI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax/
MERGED_SQANTI_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/7_sqanti/sqanti_relax_merged
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/target_gene_annotation

##-------------------------------------------------------------------------

chromosomes=($(printf "chr%s " $(seq 1 22) X Y))

cd ${MERGED_SQANTI_DIR}
for chrNum in ${chromosomes[@]}; do 
	echo $chrNum
	if compgen -G "${SQANTI_DIR}/*${chrNum}*RulesFilter_result_classification.txt" > /dev/null; then

		if [ -f WholeTargeted_RulesFilter_${chrNum}_result_classification_filenames.txt ]; then

			echo "$chrNum already merged"

		else 
			# merge classification file
			ls ${SQANTI_DIR}/*${chrNum}*RulesFilter_result_classification.txt > WholeTargeted_RulesFilter_${chrNum}_result_classification_filenames.txt

			# FNR==1 && NR!=1 { next }: This skips the first row (FNR==1) for all files except the first (NR!=1) file
			awk 'FNR==1 && NR!=1 { next } { print }' $(cat WholeTargeted_RulesFilter_${chrNum}_result_classification_filenames.txt) > WholeTargeted_RulesFilter_${chrNum}_result_classification.txt
			awk '$NF == "Isoform"' WholeTargeted_RulesFilter_${chrNum}_result_classification.txt > WholeTargeted_RulesFilter_${chrNum}_result_classification_isoform.txt

			# merge fasta
			ls ${SQANTI_DIR}/*${chrNum}*_corrected.fasta > WholeTargeted_${chrNum}_correctedfasta_filenames.txt
			cat $(cat WholeTargeted_${chrNum}_correctedfasta_filenames.txt) > WholeTargeted_${chrNum}_corrected.fasta

			# merge gtf
			ls ${SQANTI_DIR}/*${chrNum}*corrected.gtf > WholeTargeted_${chrNum}_correctedgtf_filenames.txt
			cat $(cat WholeTargeted_${chrNum}_correctedgtf_filenames.txt) > WholeTargeted_${chrNum}_corrected.gtf

			# merge gff
			ls ${SQANTI_DIR}/*${chrNum}*corrected.gtf.cds.gff > WholeTargeted_${chrNum}_correctedgff_filenames.txt
			cat $(cat WholeTargeted_${chrNum}_correctedgff_filenames.txt) > WholeTargeted_${chrNum}_corrected.gtf.cds.gff

		fi

	else

		echo $chrNum still collapsing

	fi

done
