## --------------------------------------
#!/bin/bash

# paths
pbmm2Aligned=/media/disk2/sfari_RB/5_cupcake/5_align/
DEMUX_DIR=/media/disk2/sfari_RB/7_demux/
MERGED_CHROM_DIR=/media/disk2/sfari_RB/5_cupcake/mergedChrom

export PATH=$PATH:~/installs/LOGen/assist_ont_processing/

mamba activate SQANTI3.env

# set batch array (run per chromosome)
chromosomes=($(printf "chr%s " $(seq 1 22) X Y))
for i in {0..23}; do 
	echo $i	
	chrNum=${chromosomes[${i}]}  
	echo $chrNum
	demux_cupcake_collapse.py ${MERGED_CHROM_DIR}/${chrNum}.read_stat.renamed.txt ${DEMUX_DIR}/WholeTargeted_sample_id.csv -o WholeTargeted_demux_${chrNum} --dataset ont -d ${DEMUX_DIR}
done
