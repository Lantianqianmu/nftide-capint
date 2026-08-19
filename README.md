
# Capint Pipeline #
Nextflow pipeline for detecting HBV integration from HBV-probe enriched libraries

## System requirements ##
The pipeline relies on _Perl 5.32.1_, and therefore runs only on a UNIX-based system. We have tested the pipeline on _Ubuntu 22.04_ server in our lab.  
The pipeline requires >= 64 GB memory. It is recommended to have at least 8 cores. We typically run the pipeline on a server with dual Epyc 7542 and 512 GB memory.

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10.2
openjdk | 17.0.8
Perl | 5.32.1
python | 3.12.13
cutadapt | 5.2
bwa | 0.7.19
bowtie2 | 2.5.5
pigz | 2.8
samtools | 1.23.1
bbmap | 39.81

The environment is compatible with nftide-caphic pipeline. If the environment for nftide-caphic has been set, you can use it directly, and install bbmap.

## Installation ##
Installation should finish within 10 minutes (Epyc 7542).  
(1) Setup the environment using conda or mamba:  
```
mamba create -n capint -c bioconda -c conda-forge python=3.12.13 openjdk=17.0.8 nextflow=25.10.2 cutadapt=5.2 bzip2=1.0.8 bowtie2=2.5.5 bwa=0.7.19 pigz=2.8 samtools=1.23.1 bbmap=39.81
```    
Then activate the conda environment:  
```
mamba activate capint
```
(2) Clone the repository with `git clone`, and execute
```
cd nftide-capint
```
(3) You also need to prepare bwa index files for host genome. Please refer to their manuals to generate indexed genome files.

## Overview ##
This pipeline detects HBV integration from HBV-probe captured libraries. The pipeline is inspired by the HIVID method (PMID: 23867110). It first used `bwa mem` to align reads against the HBV and host genome. Then, reads that can be mapped to both genome (chimeric reads) are extracted, and PE-reads are assembled into single contigs with `bbmerge.sh`. Both the assembled and unassembled chimeric reads are remapped to the HBV and host genome. Finally, breakpoints are called from the remapping results.

## Usage ##
(1) Prepare the `samplesheet.csv`. The csv file __must__ contain 3 columns with defined column names:  
`sample`: Name of the sequenced library. For example, `demo-1`. It will be the prefix of the output. Note: Different fastqs with same sample name will be merged before processing.  
`fastq_1`: Path to read 1.  
`fastq_2`: Path to read 2.  

(2) Unzip the files in the genome folder. It contains genome fasta files for HBV. Build bwa mem index for the HBV genome. This will be used to identify the optimal HBV genome for bwa to align against for the capture HiC library.  

(2) Unzip the files in the genome folder. It contains genome fasta files for HBV. Build bwa mem index for the HBV genome. This will be used to identify the optimal HBV genome for bwa to align against for the capture HiC library.  

(3-1) If the genotype of HBV is unknown, run the nftide-caphic pipeline. Change directory to nftide-caphic with `cd nftide-caphic`, and execute:
```
nextflow run nftide-caphic.nf \
  -output-dir your_output_dir \
  --input_csv samplesheet.csv \
  --hostGenome hg38 \
  --original_chromsizes hg38_original_chrom_sizes.tsv \
  --scaled_chromsizes hg38_scaled_chrom_sizes.tsv \
  --host_fa path_to_host_fasta \
  --allHBV_fa path_to_HBV_fasta \
  --allHBV_bwaindex path_to_HBV_bwa_index \
  --bin_size 1000 \
  --covscore true \
  --readcounts true \
  --calc_accurate_coord false \
  -with-report nftide-caphic_report.html \
  -with-timeline nftide-caphic_timeline.html \
  -bg -resume
```
`-output-dir`: Path to the output directory.  
`--input_csv`: Path to samplesheet.csv as described in **step (1)**.  
`--hostGenome`: the `--assembly` parameter for `pairtools parse`.  
`--original_chromsizes`: Path to `hg38_original_chrom_sizes.tsv`, which contains original chromsizes of the genome.  
`--scaled_chromsizes`: Path to `hg38_original_chrom_sizes.tsv`, which contains scaled chromsizes of the genome.  
`--host_fa`: fasta file of the host genome. Used to create merged host-virus bwa index.  
`--allHBV_fa`: fasta file of the HBV genome, which is acquired from **step(2)**. Used to create merged host-virus bwa index.  
`--allHBV_bwaindex`: bwa index of the HBV genome, which is acquired from **step(2)**. Used to decide the optimal HBV genome.  
`--bin_size`: Used in `cooler cload pairs`.  
`--covscore`: If true, the HBV genome with maximum mapped_coverage * mapped_read_counts will be selected as the optimal genome. Note: if mapped_read_counts > 100 & mapped_coverage < 20, the genome will not be considered since there might be biased amplification (high reads but low coverage). Default: true. Valid options: true, false.  
`--readcounts`: If true, the HBV genome with maximum mapped_read_counts will be selected as the optimal genome. Default: true. Valid options: true, false. Compatible with `--covscore` as a new instance will be initiated in parallel.  
`--calc_accurate_coord`: Whether mapping coordinates for HiC pairs should be calculated when generating contact beds and matrices. Default: false. Valid options: true, false.   


## Expected output ##
Go to `-output-dir`. The pipeline will create folders named according to the `sample` column in you csv file. Each folder contain 7 subfolders:  
`fastqs`: Merged and adapter-trimmed fastqs.  
`meta`: Optimal and assembled HBV genomes, and qc metrices.  
`bams`: Aligned and filtered bam files.  
`pairs`: pairtools outputs.  
`beds`: tsv files containing contact pairs. The xbed contains contact pairs in the following order: `chr, start, end`.  
`cools`: cool and mcool files.  
`contact_mats`: 100kb- and 1mb- binned contact matrices in dataframe format, acquired from `cooler dump`.  
`bw`: virus-host contact and virus-virus contact bw files for genome browser.  












