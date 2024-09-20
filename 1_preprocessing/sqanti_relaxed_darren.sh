## --------------------------------------
#!/bin/bash

# paths
UTILS_DIR=/home/darren/installs/SFARI_developmentalgenomics/1_preprocessing
MERGED_CHROM_DIR=/media/disk2/sfari_RB/5_cupcake/mergedChrom
SPLIT=/media/disk2/sfari_RB/6_sqanti/split
SQANTI_OUTPUT=/media/disk2/sfari_RB/6_sqanti/sqanti_relaxed

## Software 
export SOFTDIR=/home/darren/installs/
export CUPCAKE=${SOFTDIR}/cDNA_Cupcake
export SEQUENCE=$CUPCAKE/sequence
export SQANTI3_DIR=${SOFTDIR}/SQANTI3-5.2.2
export CAGE_PEAK=${SOFTDIR}/SQANTI3-5.2.2/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=${SOFTDIR}/SQANTI3-5.2.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt
export refAnno=/media/disk2/sfari_RB/references/gencode.v38.annotation.gtf
export refFile=/media/disk2/sfari_RB/references/hg38.fa
export SQANTI_JSON=/home/darren/installs/SFARI_developmentalgenomics/1_preprocessing/sqanti_filter_adapted.json

source activate SQANTI3.env
##-------------------------------------------------------------------------

cd ${SQANTI_OUTPUT}

## 335 lines in ${SPLIT}/WholeTargeted_chunks.txt
for num in {1..355}; do
	echo $num
	# print the SLURM_ARRAY_TASK_ID nth line and delete the rest
	i=`sed -n "${num}p;d" ${SPLIT}/WholeTargeted_chunks.txt`
	prefix=$(basename $i .txt)

	## extract variables from $i for downstream labelling 
	## e.g. i=WholeTargeted_chr10_chunk00.txt
	## chr="chr10_chunk00.txt" (extract after WholeTargeted)
	chr=${i#WholeTargeted_}

	## chr_chunk=chunk10_chunk00 (remove txt before chr10_chunk00.txt)
	chr_chunk=${chr%.txt}

	## remove everything after (and including) the first underscore _ in the string of chr10_chunk00.txt
	## chr=chr10
	chr=${chr%%_*}

	echo "Processing $i; output chromosome = $chr; output gff: $prefix.gff; output sqanti: $chr_chunk"
	## extract the relevant gff to run sqanti
	## renamed.gff from renaming the PB to ONT<chr> in the gff after cupcake collapse
	grep -wF -f ${SPLIT}/${i} ${MERGED_CHROM_DIR}/${chr}.renamed.gff >> ${SPLIT}/${prefix}.gff

	## run sqanti
	python ${SQANTI3_DIR}/sqanti3_qc.py ${SPLIT}/${prefix}.gff ${refAnno} ${refFile} -o WholeTargeted_collapsed${chr_chunk} \
	  --report skip --genename --skipORF \
	  --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA}
	  
	## run sqanti filter
	python ${SQANTI3_DIR}/sqanti3_filter.py rules \
	  --output WholeTargeted_collapsed_filtered${chr_chunk} \
	  --skip_report --gtf WholeTargeted_collapsed${chr_chunk}_corrected.gtf \
	  --json_filter ${SQANTI_JSON} WholeTargeted_collapsed${chr_chunk}_classification.txt
done
