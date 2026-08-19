
# Capint Pipeline #
Nextflow pipeline for detecting HBV integration from HBV-probe enriched libraries

## System requirements ##
The pipeline relies on _Perl 5.32.1_, and therefore runs only on a UNIX-based system. We have tested the pipeline on _Ubuntu 22.04_ server in our lab.  
The pipeline requires >= 32 GB memory. It is recommended to have at least 8 cores. We typically run the pipeline on a server with dual Epyc 7542 and 512 GB memory.

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

(2) Prepare the per-sample virus genome `meta_all_assembled_fa.csv`. The csv file __must__ contain 2 columns with defined column names:  
`sample`: Name of the sequenced library. Must be the same with the sample name in step (1).
`sequence`: Genome sequences.   

(3) Change directory to nftide-caphic with `cd nftide-capint`, and execute:
```
nextflow run main.nf \
  -output-dir your_output_dir \
  --input_csv samplesheet.csv \
  --HBVfa_csv meta_all_assembled_fa.csv \
  --host_bwa_dir path_to_host_bwa_dir \
  --host_bwa_prefix host_bwa_prefix \
  --bwaMinScore 15 \
  --minMatch 30 \
  --minOverlap 6 \
  --maxMismatchRate 0.2 \
  --mergeDist 20 \
  --peAssembleCustom false \
  --keepUnmerged true \
  --supportDist 300 \
  -with-report nftide-caphic_report.html \
  -with-timeline nftide-caphic_timeline.html \
  -bg -resume
```
`-output-dir`: Path to the output directory.  
`--input_csv`: Path to samplesheet.csv as described in **step (1)**.  
`--HBVfa_csv`: Path to meta_all_assembled_fa.csv as described in **step (2)**.  
`--host_bwa_dir`: Directory to bwa index of the host genome.  
`--host_bwa_prefix`: Prefix of the bwa index of the host genome.  
`--bwaMinScore`: Value for `bwa mem -T`. Minimum score to report an alignment (0 = too slow/huge output). Default: 15.  
`--minMatch`: Minimum match length (bp) on host & HBV for a breakpoint. Default: 30. Do not set this value lower that 25, as it will introduce huge amount of false positive due to homologous sequences between HBV and host.  
`--minOverlap`: Minimum overlap required for pair-end assembly. Default: 6.  
`--maxMismatchRate`: Maximum mismatch rate allowed in the overlap. Default: 0.2.  
`--mergeDist`: Merge breakpoints within this interval. Default: 20.  
`--peAssembleCustom`: true = use the custom PE assembly script, false (default) = use `bbmerge.sh`.  
`--keepUnmerged`: true (default) = parse and count non-overlapping chimeric pairs when calling breakpoints; false = ignore them.  
`--supportDist`: Maximum distance (bp) from an unmerged support pair to a breakpoint (host side). Default: 300.  


## Expected output ##
Go to `-output-dir`. The pipeline will create folders named according to the `sample` column in you csv file. Each folder contain 7 subfolders:  
`fastqs`: Merged, adapter-trimmed, and chimeric-reads filtered fastqs.  
`qc`: Proportion of virus reads for each library.  
`assembly`: PE-assembled and unassembled reads.  
`sam`: SAM files that generated from remapping chimeric reads to the host and virus genome.  
`breakpoints`: Breakpoints called.  












