#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --mem=200G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/
QCRoot=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC
cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC
mkdir -p 1b_demultiplex_merged 3_minimap 5_cupcake
mkdir -p 5_cupcake/5_align/PAF 5_cupcake/6_collapse 5_cupcake/7_sqanti3

module load Miniconda2
source activate lrp 

# stats for demuxed files (raw file) 
demuxDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/1_raw/PromethION/Project_10380/Upload_to_SRA
demuxOutDir=${QCRoot}/1b_demultiplex_merged
for file in ${demuxDir}/*fastq*; do 
	gval=$(basename $file .fastq.gz)
	echo "Processing ${file}, output: $gval"
  if [ -s ${demuxOutDir}/${gval}_readstats.txt ]; then
    echo "Already processed"
  else
    echo "Running"
  	#seqkit stats -a ${file} > ${demuxOutDir}/${gval}_readstats.txt
   fi 
done

# stats for aligned files (3_minimap2)
alignedDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/3_aligned
alignedOutDir=${QCRoot}/3_minimap/original
mkdir ${alignedOutDir} 
for file in ${alignedDir}/*sorted.sam*; do 
	name=$(basename $file .fastq.sorted.sam)
	echo "Processing ${file}, output: $name"
	htsbox samview -pS $file > $alignedOutDir/${name}.paf
	awk -F'\t' '{if ($6!="*") {print $0}}' $alignedOutDir/${name}.paf > $alignedOutDir/${name}.filtered.paf
	awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $alignedOutDir/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $alignedOutDir/${name}"_mappedstats.txt"
done
cd /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/1_preprocessing/1_QC
python $LOGEN_ROOT/miscellaneous/replace_filenames_with_csv.py \
--input=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC/3_minimap/original \
--dir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC/3_minimap/ --copy --file rename_sfari_barcode.csv --ext=txt

# stats for further aligned files 5_align
alignedDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/A_Whole/5_isoseq/pbmm2_align
alignedOutDir=${QCRoot}/5_cupcake/5_align/PAF
for file in ${alignedDir}/*_aligned_clean_aligned.bam; do 
	name=$(basename $file _aligned_clean_aligned.bam)
	echo "Processing ${file}, output: $name"
	htsbox samview -pS $file > $alignedOutDir/${name}.paf
	awk -F'\t' '{if ($6!="*") {print $0}}' $alignedOutDir/${name}.paf > $alignedOutDir/${name}.filtered.paf
	awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $alignedOutDir/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $alignedOutDir/${name}"_mappedstats.txt"
done
cd /lustre/projects/Research_Project-MRC148213/lsl693/scripts/SFARI_developmentalgenomics/1_preprocessing/1_QC
python $LOGEN_ROOT/miscellaneous/replace_filenames_with_csv.py \
--input=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC/5_cupcake/5_align/PAF/original \
--dir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/QC/5_cupcake/5_align/PAF/ --copy --file rename_sfari_id.csv --ext=txt


# sqanti output
sqantiDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/C_Whole_Targeted/9_sqanti_final
sqantiOutDir=${QCRoot}/5_cupcake/7_sqanti3

cp ${sqantiDir}/*RulesFilter* ${sqantiOutDir}
cp ${sqantiDir}/*filtering_reasons* ${sqantiOutDir}


