# nftide-capint Pipeline #

Nextflow pipeline for HBV integration detection from paired-end HBV-probe capture libraries.

## System requirements ##

The pipeline runs on Linux with Bash, Python 3, Nextflow and Java. Bioinformatics
executables must be available on the task PATH. This workflow does not activate
a separate Conda environment automatically; use an environment containing the
software below or explicitly configure PATH before launching Nextflow.

Plan storage for the host BWA index, per-sample HBV indexes, intermediate FASTQs
and alignments, and published copies of results. A fixed minimum memory
requirement has not been benchmarked. The current `nextflow.config` requests
16 CPUs per process and permits up to six concurrent tasks per process.
Reduce concurrency on smaller servers. Published results use copy mode;
`cleanup = false` retains work files for resuming and diagnosis.

## Software dependencies ##

Dependencies | Version
------------- | -------------
Nextflow | 25.10.2
Java / OpenJDK | 23.0.2
Python | 3.12.13
BWA | 0.7.19
samtools | 1.23.1
cutadapt | 5.2
BBMerge / BBMap | 39.81

Dependencies are resolved from the active environment. The Python analysis
scripts use the standard library and are provided in `bin/`.

## Installation ##

(1) Enter the pipeline directory:

```bash
cd /data/xrz/capint/nextflow
```

(2) Create and activate a Conda environment with explicit package versions:

```bash
conda create -n nftide-capint \
  --override-channels -c conda-forge -c bioconda \
  python=3.12.13 nextflow=25.10.2 openjdk=23.0.2 \
  bwa=0.7.19 samtools=1.23.1 cutadapt=5.2 bbmap=39.81
conda activate nftide-capint
```

`bbmap` supplies `bbmerge.sh`; use the `bwa` package for BWA-MEM. These versions
come from the existing local installations. The combined environment creation
command has not been freshly solved or installed as part of this documentation
update. Conda will resolve the remaining dependencies; this is not a complete
build-level lock file.

Check the executables before launching:

```bash
python --version
java -version
nextflow -version
samtools --version
cutadapt --version
command -v bwa
command -v bbmerge.sh
```

Keep this environment active when running `nextflow run main.nf`; processing
tasks inherit PATH. The pipeline does not provide a pinned environment file
or automatically activate an environment. After installation, record the exact
resolved packages for reproducibility:

```bash
conda list -n nftide-capint --explicit > nftide-capint-conda-explicit.txt
```

This repository runs the HIVID-derived BWA/BBMerge workflow described below.

(3) Prepare a human BWA index and configure its directory and prefix. Current
defaults are:

```text
Human BWA index directory: /data/xrz/ref/hg38/hg38_bwa
Human BWA index prefix:    hg38.fa
HBV reference CSV:         /data/xrz/capint/nextflow/meta_all_assembled_fa.csv
```

`host_bwa_prefix` is the index basename within `host_bwa_dir`, not the directory
itself. The pipeline extracts each sample's HBV sequence from the CSV and
builds its HBV BWA index automatically.

## Overview ##

The pipeline maps with BWA-MEM, reconstructs overlapping pairs, calls sequenced
human–HBV junctions, and adds compatible discordant-pair evidence.

FASTQs with the same sample name are merged before adapter trimming. Candidate
pairs are selected from separate human and sample-specific HBV alignments,
sequence-deduplicated, and assembled with BBMerge by default. The optional
Python assembler can be selected with `peAssembleCustom`.

Sequenced junctions take priority when assigning discordant-pair support.
Pairs without a compatible sequenced junction retain inward-facing mate-end
coordinates, explicitly marked as estimates.

## Usage ##

(1) Prepare the input CSV and per-sample HBV sequence CSV.

The samplesheet must contain the following columns:

`sample`: Sample name and output folder; rows sharing a sample name are merged.  
`fastq_1`: Path to the gzipped read 1 FASTQ.  
`fastq_2`: Path to the matching gzipped read 2 FASTQ.

```csv
sample,fastq_1,fastq_2
sample1,/path/sample1_lane1_R1.fastq.gz,/path/sample1_lane1_R2.fastq.gz
sample1,/path/sample1_lane2_R1.fastq.gz,/path/sample1_lane2_R2.fastq.gz
sample2,/path/sample2_R1.fastq.gz,/path/sample2_R2.fastq.gz
```

The HBV CSV must contain `sample` and `sequence` columns, with one complete HBV
sequence per sample. Every input sample must have an HBV reference entry. `sequence` contains the
nucleotide sequence itself, not a FASTA filename. For example:

```csv
sample,sequence
sample1,ACGTACGTACGTACGT
sample2,ACGTACGTACGTACGT
```

These short sequences illustrate the format only; supply the actual assembled
HBV genomes. Use simple sample names without spaces or shell-special characters.
HBV positions refer to each sample's supplied sequence.

(2) Run `main.nf` with the input CSV, per-sample HBV sequence CSV and prebuilt
human BWA index:

```bash
nextflow run main.nf \
  --input_csv samplesheet.csv \
  --HBVfa_csv meta_all_assembled_fa.csv \
  --host_bwa_dir /data/xrz/ref/hg38/hg38_bwa \
  --host_bwa_prefix hg38.fa \
  -output-dir /data/xrz/capint/output_test \
  -with-report nftide-capint_report.html \
  -with-timeline nftide-capint_timeline.html \
  -with-trace nftide-capint_trace.tsv \
  -resume -bg
```

### General parameters ###

Parameter | Default in the current pipeline | Description
------------- | ------------- | -------------
`-output-dir` | `/data/xrz/capint/output_test` | Published results directory; each sample receives its own subdirectory. This is a Nextflow option, not `--output_dir`.
`--input_csv` | `/data/xrz/capint/nextflow/samplesheet.csv` | FASTQ samplesheet with `sample`, `fastq_1`, and `fastq_2` columns.
`--HBVfa_csv` | `/data/xrz/capint/nextflow/meta_all_assembled_fa.csv` | Sample-specific HBV sequences in `sample` and `sequence` columns.
`--host_bwa_dir` | `/data/xrz/ref/hg38/hg38_bwa` | Directory containing the prebuilt host BWA index.
`--host_bwa_prefix` | `hg38.fa` | Index basename within that directory, including the FASTA extension when it is part of the index name. For example, sidecars are named `hg38.fa.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
`-resume` | Disabled unless supplied | Reuse eligible completed tasks; retain both `work/` and `.nextflow/`.
`-bg` | Disabled unless supplied | Run Nextflow in the background; inspect `.nextflow.log` for progress.
`-with-report` | Disabled unless supplied | Write an HTML execution report to the specified path.
`-with-timeline` | Disabled unless supplied | Write an HTML task timeline to the specified path.
`-with-trace` | Disabled unless supplied | Write a tab-delimited task execution trace to the specified path.

Nextflow options use one hyphen; pipeline parameters use two. Resource settings
are configured in `nextflow.config`, not through a custom `--cpus` parameter.

### Alignment and assembly parameters ###

Analysis parameters are defined in `main.nf` and can be overridden with
`--parameter value`. Output and process-resource defaults are in
`nextflow.config`. Each process requests 16 CPUs, with `maxForks=6`; adjust
these values for available resources and concurrent stages.

Parameter | Default | Description
------------- | ------------- | -------------
`bwaMinScore` | `15` | BWA-MEM alignment-score reporting threshold (`-T`)
`minOverlap` | `6` | Minimum paired-end assembly overlap in bp
`maxMismatchRate` | `0.2` | Overlap mismatch setting passed to the selected assembler
`peAssembleCustom` | `false` | Use the optional Python assembler instead of BBMerge
`keepUnmerged` | `true` | Retain unmerged reads for junction calling and discordant-pair support
`nnssThreshold` | `1` | Deprecated; the final table is not filtered by this value

### Junction and support parameters ###

The BWA-MEM junction and clustering parameters are:

Parameter | Default | Description
------------- | ------------- | -------------
`minMatch` | 30 | Each anchor must align strictly more than 30 query bases by default; deletions do not count
`minUnique` | 20 | Minimum aligned bases exclusive to each reference
`maxJunctionOverlap` | 10 | Maximum human/HBV query overlap
`junctionScoreMargin` | 10 | Score advantage for overlapping junction anchors; also used to check competing discordant-pair mappings
`minMapq` | 20 | Threshold for uniquely mapped support
`altScoreDelta` | 5 | Score window for alternative junction explanations
`mergeDist` | 20 | Maximum clustering distance on each reference
`supportDist` | 500 | Maximum total human + HBV unsequenced span for pair assignment

`supportDist` limits assignment of an unsequenced pair interval to a sequenced
junction. The host and HBV distances are measured towards that interval using
the reported strands. Each must lie between -10 and 500 bp, and their sum must
be at most 500 bp at the current default. It also clusters unresolved estimated
pairs within 500 bp on **each** reference; this is distinct from `mergeDist=20`
for sequenced calls. The standalone `count_support.py` default remains 300;
Nextflow explicitly passes 500. Increasing this setting does not determine an
exact junction position inside an unsequenced interval.

Discordant anchors must cover at least 90% of their reads and lack a comparable
opposite-genome explanation. Pair assignment allows up to 10 bp of boundary
slop. These are initial analysis settings, not empirically calibrated false
positive thresholds. No merging through shared repeat alternatives is performed.

There is **no maximum sequenced-gap length**. Disjoint human and HBV anchors
are accepted when each has more than `minMatch` aligned query bases (strictly
more than 30 bp by default), regardless of the intervening sequence length.
Overlapping alignments still require the overlap, exclusive-sequence and score
checks; alignment ambiguity remains annotated.

The HBV reference remains linear and sample-specific; this revision does not
standardize HBV coordinates between samples or resolve circular-origin crossings.

### Validation and rerunning ###

Run `python -m unittest discover -s tests -v` for regression tests.
With all dependencies on PATH, `python tests/smoke_pipeline.py` runs the actual
Nextflow workflow against deterministic tiny synthetic references in `/tmp`.
Pass `--custom` to exercise the optional Python assembler instead of BBMerge.

The breakpoint schema and support semantics changed. Regenerate candidate
FASTQs and every downstream stage; old breakpoint tables cannot be mixed with
new ones. Prefer a fresh output directory for comparison with previous results.
The new BWA commands and extraction inputs also invalidate their old task cache
entries. Do not interpret older results as corrected without rerunning.

### BBMerge wrapper output contract ###

`assemble_pairs.py` is the default `PE_ASSEMBLE` wrapper and runs BBMerge itself.
It produces the same consensus sequences and qualities as the BBMerge run it
invokes. The final FASTQs normalize only read names: merged reads have no mate
suffix; unmerged mates have `/1` and `/2`. Molecule IDs and copy counts remain
in the names. Gzip bytes and record order need not match a separate invocation.
The native merged output is retained in the task directory as
`<sample>_bbmerge_merged.fq.gz` for comparison.

After checking read conservation, it stages the final merged and two unmerged
FASTQs and moves them into their expected filenames. `keepUnmerged=false`
produces valid empty gzip files for the two unmerged outputs, while the QC
still reports the actual unmerged count and `emitted_unmerged_pairs=0`.

This wrapper is distinct from `pe_assemble.py`, the optional custom Python
assembler selected by `peAssembleCustom=true`. That custom algorithm is not
expected to reproduce BBMerge's consensus decisions or qualities.


### Custom assembler (`peAssembleCustom=true`) ###

`pe_assemble.py` emits the same merged FASTQ, separate unmerged R1/R2 FASTQs,
assembly log and assembly QC JSON as the BBMerge branch. Molecule/copy tags
are preserved; merged names have no mate suffix and unmerged names have
`/1` and `/2`. Input pairs are streamed in matching order; missing or
mismatched mates cause an error rather than silent read loss. With
`keepUnmerged=false`, both unmerged files are valid empty gzip files, while
QC still reports the number of unmerged pairs.

The custom algorithm reverse-complements R2 and chooses the longest ungapped
R1-suffix/R2-prefix overlap satisfying `minOverlap` and `maxMismatchRate`.
Ambiguous bases count as mismatches. Unmerged reads and nonoverlapping tails
retain their original qualities. Matching overlap bases take the higher
observed quality; disagreements take the higher-quality base and the absolute
Phred-quality difference, with equal-quality disagreements becoming `N` at
Phred 0. These are simple consensus rules, not BBMerge's quality model.


## Expected output ##

The pipeline creates a separate folder for each sample under `-output-dir`:

```text
output/<sample>/
├── fastqs/
│   ├── <sample>_cutadapt_R1.fq.gz
│   ├── <sample>_cutadapt_R2.fq.gz
│   ├── <sample>_cutadapt.log
│   ├── <sample>_first_align.log
│   ├── <sample>_chimeric_R1.fq.gz
│   ├── <sample>_chimeric_R2.fq.gz
│   ├── <sample>_chimeric.log
│   └── <sample>_hbv.fa
├── assembly/
│   ├── <sample>_assembly.log
│   ├── <sample>_assembly_qc.json
│   ├── <sample>_merged.fq.gz
│   ├── <sample>_unmerged_R1.fq.gz
│   └── <sample>_unmerged_R2.fq.gz
├── sam/
│   ├── <sample>_host.bam
│   ├── <sample>_hbv.bam
│   ├── <sample>_remap_host.sam
│   ├── <sample>_remap_hbv.sam
│   ├── <sample>_support_host.sam
│   └── <sample>_support_hbv.sam
├── qc/
│   ├── <sample>_hbv_qc.tsv
│   └── <sample>_molecule_map.tsv
└── breakpoints/
    ├── <sample>_breakpoints.txt
    ├── <sample>_breakpoints_combined.txt
    ├── <sample>_breakpoints_merged.txt
    ├── <sample>_breakpoints_final.txt
    ├── <sample>_breakpoints_primary.txt
    └── <sample>_nnss_summary.txt
```

### FASTQs and alignments ###

File suffix | Description
------------- | -------------
`_cutadapt_R1.fq.gz`, `_cutadapt_R2.fq.gz` | Adapter-trimmed paired FASTQs after merging all lanes for the sample. These supply original sequences/qualities and the NNSS denominator.
`_cutadapt.log` | Adapter trimming and retained/discarded read statistics. Fixed settings use adapter `CTGTCTCTTATACACATCT` on both mates, minimum length `2:2`, overlap 1, and pair filtering `any`.
`_first_align.log` | Log of the initial separate host/HBV mappings.
`_chimeric_R1.fq.gz`, `_chimeric_R2.fq.gz` | Sequence-deduplicated candidate pairs with mapping evidence to both references. Candidate status alone does not establish a junction.
`_chimeric.log` | Candidate-pair and sequence-deduplicated molecule counts.
`_hbv.fa` | Sample-specific HBV reference used for mapping and HBV coordinates.
`_merged.fq.gz` | One consensus record per successfully assembled candidate pair; names preserve molecule/copy tags without mate suffixes.
`_unmerged_R1.fq.gz`, `_unmerged_R2.fq.gz` | Unassembled candidate pairs with `/1` and `/2` suffixes, retained for sequenced-mate calls and discordant support. Valid empty gzip files when none remain or `keepUnmerged=false`.
`_assembly.log` | BBMerge diagnostics or custom assembler summary, depending on the selected branch.
`_assembly_qc.json` | Assembly accounting, described below; emitted by both branches.
`_host.bam`, `_hbv.bam` | Initial mappings of trimmed paired reads against the separate references. These are not coordinate-sorted/indexed BAM deliverables.
`_remap_host.sam`, `_remap_hbv.sam` | Mappings of assembled consensus reads used to call sequenced junctions.
`_support_host.sam`, `_support_hbv.sam` | Mappings of retained unmerged pairs used to find within-mate junctions or discordant support.

### Capture, molecule and assembly QC ###

`qc/<sample>_hbv_qc.tsv` has a header and one data row:

Column | Meaning
------------- | -------------
`sample` | Sample/library identifier.
`total_read_pairs` | Number of trimmed R1 records, representing paired input after trimming.
`hbv_read_pairs` | Number of distinct query names with at least one reported HBV alignment; not specifically integration-supporting pairs.
`prop_hbv_reads` | `hbv_read_pairs / total_read_pairs`, printed to six decimal places; a fraction, not a percentage. Zero for empty input.

`qc/<sample>_molecule_map.tsv` has `read_id` and `molecule_id` columns.
Each original candidate pair is linked to its sequence-based representative;
several read IDs can share a molecule ID. Noncandidate pairs are not included.

`assembly/<sample>_assembly_qc.json` contains:

Key | Meaning
------------- | -------------
`input_pairs` | Deduplicated candidate pairs passed to the assembler.
`merged_pairs` | Pairs emitted as assembled consensus records.
`unmerged_pairs` | Pairs not assembled, including those excluded when `keepUnmerged=false`.
`keep_unmerged` | Boolean recording whether unassembled pairs are retained downstream.
`emitted_unmerged_pairs` | Number written to each unmerged FASTQ: `unmerged_pairs` when retained, otherwise 0.

A successful assembly accounts for `input_pairs = merged_pairs + unmerged_pairs`.
Capture QC and assembly QC describe different populations: HBV-aligned pairs
include nonchimeric pairs and duplicates, whereas assembly uses deduplicated
candidates. For N5, 14,373 HBV-aligned pairs yielded 14,314 candidate pairs and
3,749 sequence-deduplicated assembly inputs. Thus an HBV capture count alone
does not predict the merged FASTQ record count.

### Breakpoint files ###

The breakpoint files represent successive processing stages:

File suffix | Contents
------------- | -------------
`_breakpoints.txt` | Junctions detected from **merged reads**, before clustering or adding unmerged-read support.
`_breakpoints_combined.txt` | Adds junctions detected within individual unmerged mates and support from discordant pairs. Sequenced calls are clustered first. Pairs support compatible sequenced breakpoints; otherwise, they produce estimated coordinates.
`_breakpoints_merged.txt` | Clusters compatible calls using human position, HBV position and orientation. Counts each molecule proxy once per cluster.
`_breakpoints_final.txt` | Adds `efr` and `nnss` to the merged table. **No confidence filtering is applied.**
`_breakpoints_primary.txt` | Currently has the **same contents as `_final.txt`**. The filename is retained for compatibility; alternative mappings are stored within each row rather than as separate rows.

**Use `_breakpoints_final.txt` for downstream analysis.** Check the evidence
resolution, support counts and ambiguity before interpreting a call. These
file descriptions apply to outputs generated after the fixes; existing files
from earlier runs retain the old behavior.

`_nnss_summary.txt` contains one tab-separated line without a header:
`sample<TAB>efr`. It records the sample name and effective input-pair count
used as the normalization denominator. `qc/*_hbv_qc.tsv` reports total input
pairs, HBV-aligned pairs and their proportion.

### Junction evidence ###

`breakpoint_resolution` distinguishes:

- `sequenced`: adjacent human and viral alignments in a sequenced read.
- `microhomology`: limited overlap; coordinates select a consistent split within
  the overlap. The `*_ci_start/end` columns retain the coordinate ambiguity.
- `sequenced_gap`: a sequenced segment between the aligned anchors, with no
  maximum allowed length.
- `estimated_pair`: the junction is not sequenced. Coordinates are the
  inward-facing matched bases of the two mates, not exact junction positions.
  For these estimates, `*_ci_start/end` only describe the observed border range
  of clustered pairs; they are not confidence bounds for the true junction.

Unmerged reads with a junction within either mate are treated as sequenced
junction evidence. All sequenced calls are clustered first. Discordant pairs
are then assigned to a compatible sequenced call using **both references**,
orientation, and the direction towards the unsequenced interval. If no such
call exists, inward-facing mate-end coordinates are retained as estimates.
If several calls are compatible, support is assigned once to the closest and
marked ambiguous; competing event IDs remain in its molecule metadata.

Coordinates are 1-based. `host_strand` and `hbv_strand` describe traversal
**from human into HBV**, independent of read/mate orientation. This convention
allows opposite reads of the same junction to cluster consistently.

### Counting and ambiguity ###

Original trimmed FASTQs supply sequences and qualities for candidate pairs;
SAM reverse-strand sequence is never passed directly to BBMerge. Identical
full candidate-pair sequences (including swapped-mate representation)
are collapsed before assembly. Original read IDs are retained in the published
`qc/*_molecule_map.tsv`. Each representative carries a molecule ID and copy
count through assembly and remapping.

`nss` counts distinct sequence-based molecule proxies; `raw_read_pairs` retains
the pre-collapse support count. `unique_nss`, `ambiguous_nss`, `split_nss` and
`pair_nss` explain its composition. Exact-sequence deduplication is **not UMI
counting**: PCR copies with sequencing differences can remain separate, and
independent identical fragments can collapse. These counts do not establish
physical molecule independence. The workflow has no input UMI/library schema.

MAPQ below `minMapq` or competing placements yields ambiguous evidence.
Low-MAPQ host sequences are remapped with `bwa mem -a -Y`, preserving original
query orientation and CIGAR geometry. Alternative junction placements within
`altScoreDelta` are stored in `secondary_events` JSON; they are not extra
integration events. BWA's found alternatives are not guaranteed exhaustive.
Ambiguous support never becomes unique merely because it clusters with a call.

`molecules` is JSON provenance used for deduplication across stages. The final
and historical `primary` tables both retain ambiguous candidates. `event_type`
is `ambiguous` when there is no uniquely mapped supporting molecule; otherwise
it is `primary`. Neither label is a biological validation claim.

NNSS remains NSS per million trimmed input pairs (before candidate deduplication).
The final table is an unfiltered candidate table; the historical NNSS threshold
parameter remains deprecated.

### Breakpoint table columns ###

The `.txt` breakpoint files are tab-delimited tables with a header, including
when no calls are present. The raw, combined and merged stages share the evidence columns
below; `efr` and `nnss` are appended only in the final and primary tables.
JSON fields use TSV/CSV quoting: use a tab-delimited CSV reader before decoding
the JSON, rather than treating doubled quote characters as literal JSON.

Column | Meaning
------------- | -------------
`host_chr` | Host reference contig/chromosome name.
`host_pos` | Representative 1-based host anchor boundary; for `estimated_pair`, the inward-facing aligned mate end.
`host_strand` | Host traversal orientation (`+` or `-`) from human into HBV, not simply the original read's strand.
`hbv_pos` | Representative 1-based boundary on the sample-specific linear HBV reference.
`hbv_strand` | HBV traversal orientation (`+` or `-`) in the same human-to-HBV direction.
`nss` | Number of distinct sequence-based molecule proxies supporting this row.
`source` | `merged` or `unmerged` origin. A cluster is labeled `merged` when it includes merged-read evidence; this is not a complete breakdown of its support.
`host_mapq`, `hbv_mapq` | BWA MAPQ values for the representative host and HBV anchors; not aggregate confidence scores for all supporting reads.
`event_id` | Initial sample-plus-molecule identifier retained from the cluster representative; not a globally stable genomic site identifier across reruns or parameter changes.
`event_type` | `primary` when at least one supporting molecule is uniquely mapped under the current rules; otherwise `ambiguous`. Neither means experimentally validated.
`breakpoint_resolution` | `sequenced`, `microhomology`, `sequenced_gap`, or `estimated_pair`; see Junction evidence above.
`host_ci_start`, `host_ci_end` | Inclusive 1-based host boundary range retained across the evidence/cluster. Describes microhomology or observed boundary variation, not a statistical confidence interval.
`hbv_ci_start`, `hbv_ci_end` | Corresponding inclusive HBV boundary range. For estimated pairs, neither reference range bounds the unknown true junction.
`unique_nss` | Supporting molecule proxies without detected mapping ambiguity.
`ambiguous_nss` | Supporting molecule proxies with low MAPQ or competing assignments; `unique_nss + ambiguous_nss = nss`.
`split_nss` | Molecule proxies with a sequenced junction, whether in a merged read or an individual unmerged mate.
`pair_nss` | Proxies supported only by discordant mate placements; `split_nss + pair_nss = nss`. A proxy with both kinds is counted as split evidence once.
`raw_read_pairs` | Sum of copy counts for supporting proxies, before exact-sequence deduplication. This is not the total sample read count.
`secondary_events` | JSON array of alternative junction placements, including their coordinates, orientation and gap information. Empty array `[]` means none were retained, not proof of unique placement.
`molecules` | JSON object keyed by molecule ID, preserving support provenance and preventing duplicate counting across stages.
`junction_gap_length` | Sequenced query bases between the representative anchors; `0` for adjacent or overlapping anchors, `NA` for an unsequenced mate interval.
`efr` | Number of trimmed input pairs before candidate selection and sequence deduplication; normalization denominator.
`nnss` | `nss × 1,000,000 / efr`, printed to three decimal places; 0 when the denominator is 0. Reported without threshold filtering.

Each `molecules` entry records `ambiguous` (Boolean), `kind` (`split` or `pair`),
`copies` (original pair count), and `junction_gap_length`. Assigned discordant
support can additionally retain `alternative_placements` and, when several
sequenced events were compatible, `compatible_events`. The latter records
competing representative event IDs; the supporting proxy is assigned once.
For example, `nss=1` and `raw_read_pairs=2` means two identical candidate pairs
were collapsed into one sequence-based proxy, not that a read was lost.

`junction_gap_length` follows the representative breakpoint when calls are
clustered; it is not summed across reads. Individual supporting molecules retain
their own gap lengths in `molecules` JSON. Alternative placements also retain
their own gap lengths. Discordant-pair support does not replace the sequenced
representative's value. `NA` means that the pair does not sequence the gap and
its length is unknown; it must not be interpreted as zero.

The removed `maxJunctionGap` parameter and Python `--max-gap` option are no
longer used. Regenerate REMAP and downstream outputs to populate the new column.

### Large evidence fields and resuming failed runs ###

If the task reports this exact error:

```text
_csv.Error: field larger than field limit (131072)
```

Python's default CSV field limit was exceeded while reading an intermediate
breakpoint TSV. Repeat-rich candidates can contain large `secondary_events`
or `molecules` JSON fields. This error does **not** mean that the input
samplesheet is malformed. Do not shorten the JSON, remove candidate rows, or
convert the files: doing so can discard evidence.

(1) Use the corrected scripts in this repository. The fix is already included:
`bin/table_io.py` raises the CSV field limit to the largest supported value,
and `junctions.py`, `filter_signal.py`, and `collapse_secondary.py` use that
reader. If running another copy of the pipeline, update those scripts together
from this corrected copy. No extra Python package or command-line parameter
is required. Raising the limit in a separate Python terminal does not affect
the pipeline's Python processes.

(2) Keep the failed run's `work/` and `.nextflow/` directories. Return to the
**same directory from which you originally launched Nextflow**, activate the
environment, and rerun your original command with `-resume` appended. Keep the
same input files, reference paths, parameters, and output directory. For a run
originally launched from this repository, the command has this form:

```bash
conda activate nftide-capint
cd /data/xrz/capint/nextflow

# Repeat YOUR original nextflow command and add -resume:
nextflow run main.nf \
  --input_csv samplesheet.csv \
  --HBVfa_csv meta_all_assembled_fa.csv \
  --host_bwa_dir /data/xrz/ref/hg38/hg38_bwa \
  --host_bwa_prefix hg38.fa \
  -output-dir /path/to/your/original_output \
  -resume
```

Replace the example paths with those from the failed run and include any other
options you originally supplied. Nextflow reuses eligible completed tasks and
reruns failed or invalidated tasks; it does not require manual editing of
cached breakpoint tables. If you launched other runs afterwards, select the
failed run explicitly with `-resume <run-name-or-session-id>`; `nextflow log`
lists prior runs. If the work/cache directories were deleted, rerun normally;
missing cached results must be recomputed.

(3) If the **same error persists with the corrected scripts**, inspect the
new `.nextflow.log` for the failed task's work directory, then read that task's
`.command.err` and `.command.sh`. The traceback identifies the Python script
that still uses the default field limit. Check that the launch command points
to the corrected repository and that the failing reader uses `table_io.py`.
Keep the full traceback and task path for diagnosis. Other CSV errors, such as
missing columns or malformed quoting, need a separate diagnosis; increasing
the field limit does not fix them.

### Work files and reports ###

Pre-trimming lane-merged FASTQs, HBV BWA indexes, initial SAM intermediates,
remapping/support logs, and native BBMerge/interleaved temporary files remain
in Nextflow `work/` task directories; they are not all published. In particular,
`<sample>_bbmerge_merged.fq.gz` is the native BBMerge consensus output before
header normalization. Published results consume additional space because the
workflow uses copy mode.

`-with-report`, `-with-timeline` and `-with-trace` write run reports to the paths
specified in the launch command. Inspect `.nextflow.log` and the failed task's
`.command.sh`, `.command.out`, `.command.err`, and `.exitcode` for diagnostics.
Keep both `work/` and `.nextflow/` when using `-resume`.
