For more detials please go the wiki page [https://github.com/ShaodongWei/single_cell_sequencing_workflow/wiki]
# 1. Run the entire pipeline
conda config --set channel_priority flexible # set channel priority to be flexible for conda

snakemake --cores threads_number --use-conda  # conda will install dependencies into isolated environment for each step. All steps will be executed sequentially. 

snakemake --cores threads_number --use-conda --dryrun # using dryrun to check what steps will be executed


# 2. Run a specific step 
snakemake --list # to show steps

snakemake step_name --cores threads_number --use-conda # run a specific step 


# 3. Steps in the workflow 

## 3.1, demultiplexing
´´´

snakemake --cores threads_number --use-conda

´´´

This step is to demultiplex your single raw fastq files into barcoded single cell files that each file ideally represents a single droplet. You have to use 3 tandem barcode files that are allocated on the R2 file. The workflow can easily be modified to support the condition that barcode are on the R1 file or both. 

## 3.2, prune demultiplexed files based on R1 + R2 number of reads
This step is to choose the range of sequence depth, where you can specify minimal reads (leave parameter max_reads empty), maximal reads (leave parameter min_reads empty), or between minimal and maximal reads (specify both min_reads and max_reads). 

## 3.3, calculate the percentage of reads successfully demultiplexed
This step is to report the percentage of reads that are successfully demultiplexed, before applying the pruning of reads. 

## 3.4, quality control 
This step is to quality control the reads using the trimmomatic software using the default parameters for paired-end fastq files 'ILLUMINACLIP:{your_primer}:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{your.specified.trimmomatic_quality} MINLEN:30', and users have to specify your own primer used and the average per base quality score in the sliding window. 

## 3.5, mapping reads to reference 
This step is to map each barcoded fastq files to the reference. Users can choose to use bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2) as the mapper for general purposes, or the kma (https://bitbucket.org/genomicepidemiology/kma/src/master/) as the mapper when the references have similar genomes or redundant.

## 3.6, calculate depth and coverage (positions covered by reads in percentage) for each reference contig and each barcode file 
This step first calculates depth for each position in the contigs from the reference, then calculate the coverage for each reference contig. The it also reports the percentage of positions that each barcode file can cover in the reference (coverage). 

