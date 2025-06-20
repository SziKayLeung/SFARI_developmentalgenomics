### Scripts to merge whole and targeted dataset

> [!IMPORTANT]
> Given the sequencing depth of the whole dataset (n = 47 samples), had to run:
>  1. isoseq collapse by chromosome
>  2. sqanti by chromosome and chunks
>  3. obtain the transcript abundance by subsetting the isoseq collapse output file (read_stat.txt) by minimum of 2 FL reads and removing monoexonic intergenic transcripts
>  4. finally, merge sqanti output across chromosomes with abundance file

***
- 1_prepare_cupcake_wholeTargeted
    - merge the aligned bam files from whole and targeted, and split by chromosome (1 - 22, X and Y);  creates WholeTargeted.bam and WholeTargeted.sorted.bam

- 2_run_isoseqCollapse_WholeTargeted.sh
    - split bam file by chromosome and run isoseq collapse per chromosome 
    
- 3_run_sqanti_per_chromosome.sh
    - split isoseq collapsed file by 100000 lines to run sqanti in chunks

- 4_identify_nonMonoExonic_transcripts.R
    - remove monoexonic intergenic transcripts etc

- 5_filter_monoExonicTranscripts_readstat.sh
    - further filter read_stat file to remove monoexonic isoforms

- 6_retain_min2FLReadsTranscripts_readstat.sh
    - filter read_stat file minimum 2 reads in isoseq collapse readstat file —> create chr1.read_stat.renamed.min2FL.txt

- 7_demux_prep.sh
    - generate sample_id.csv for targeted and whole dataset

- 8_demux_per_sample_per_chromosome.sh
    - demultiplex per sample per chromosome using the readstat file that is now filtered for minimum 2 reads 2 samples, and removed mono-intergenic reads
      
- 9_merge_sqanti_across_chromosome.sh
    - merge sqanti output files across multiple chromosomes
      
- 10_filter_sqanti_byMonoExonic2FL2Samples.R
    - merge final sqanti output file with abundance from above (2reads, 2samples, monoexonic filtered)

> [!NOTE]
> SQANTI was generated using the relaxed json file
