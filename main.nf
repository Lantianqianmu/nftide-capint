#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// mkdir -p ~/.nextflow/lsp/v25.10
// wget https://github.com/nextflow-io/language-server/releases/download/v25.10.2/language-server-all.jar -O ~/.nextflow/lsp/v25.10/v25.10.2.jar

// nextflow run main.nf -resume -bg


params.input_csv = '/data/xrz/capint/nextflow/samplesheet.csv'


params.HBVfa_csv = "/data/xrz/capint/nextflow/meta_all_assembled_fa.csv"

// prebuilt host bwa index (shared by all samples, no need to build again)
params.host_bwa_dir = "/data/xrz/ref/hg38/hg38_bwa"
params.host_bwa_prefix = "hg38.fa"


// --- HIVID pipeline parameters ---
params.bwaMinScore = 15         // bwa mem -T: minimum score to report an alignment (0 = too slow/huge output)
params.minMatch = 30           // each host/HBV anchor must align strictly more than this many bp
params.minOverlap = 6          // PE-assembly: minimum overlap (>5 bp) required for splicing
params.maxMismatchRate = 0.2   // PE-assembly: max mismatch rate allowed in the overlap
params.mergeDist = 20          // merge breakpoints within 20 bp
params.nnssThreshold = 1       // NNSS >= 1 is regarded as a true signal
params.peAssembleCustom = false // true = use the custom PE assembly script, false (default) = use bbmerge.sh
params.keepUnmerged = true       // true (default) = count non-overlapping chimeric pairs into NSS of supported breakpoints; false = ignore them
params.supportDist = 500        // maximum total unsequenced span across host and HBV
params.minUnique = 20           // aligned bases exclusive to each reference
params.maxJunctionOverlap = 10  // maximum microhomology length
params.junctionScoreMargin = 10 // split score advantage over best single alignment
params.minMapq = 20             // lower values are retained as ambiguous evidence
params.altScoreDelta = 5        // retain competing junction placements within this score


process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."
    
    input:
    tuple val(meta), path(r1s) , path(r2s)
    
    output:
    tuple val(meta), path("${meta.id}_merged_R1.fq.gz"), path("${meta.id}_merged_R2.fq.gz"), emit: merged_fq
    
    script:
    """
    cat ${r1s.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${r2s.join(' ')} > ${meta.id}_merged_R2.fq.gz
    """
}


process CUTADAPT {
    tag "cutadapt on ${meta.id}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_cutadapt_R1.fq.gz"), path("${meta.id}_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.id}_cutadapt.log"), emit: cutadapt_log

    script:
    """
    cutadapt \
        -j ${task.cpus} -m 2:2 \
        -a "CTGTCTCTTATACACATCT" \
        -A "CTGTCTCTTATACACATCT" \
        --pair-filter=any \
        --overlap 1 \
        -o ${meta.id}_cutadapt_R1.fq.gz \
        -p ${meta.id}_cutadapt_R2.fq.gz \
        ${r1} \
        ${r2} > "${meta.id}_cutadapt.log"

    """
}


process CHECK_SAMPLES {
    tag "Checking that every sample has an HBV reference genome"

    input:
    val missing

    script:
    """
    if [ -n "${missing}" ]; then
        echo "ERROR: no HBV reference genome found for sample(s): ${missing}" >&2
        exit 1
    fi
    echo "all samples have an HBV reference genome"
    """
}


process EXTRACT_HBV_FA {
    tag "Extracting HBV genome of ${meta.id}"

    input:
    tuple val(meta), val(hbv_seq)

    output:
    tuple val(meta), path("${meta.id}_HBV.fa"), emit: hbv_fa

    script:
    """
    printf '>HBV\\n${hbv_seq}\\n' > ${meta.id}_HBV.fa
    """
}


process BUILD_HBV_INDEX {
    tag "Building bwa index for HBV genome of ${meta.id}"

    input:
    tuple val(meta), path(hbv_fa)

    output:
    tuple val(meta), path("${meta.id}_hbv_index"), emit: hbv_index

    script:
    """
    mkdir -p ${meta.id}_hbv_index
    cp ${hbv_fa} ${meta.id}_hbv_index/genome.fa
    cd ${meta.id}_hbv_index && bwa index genome.fa
    """
}


process FIRST_ALIGN {
    tag "First-time alignment (bwa mem) of ${meta.id}"

    input:
    tuple val(meta), path(r1), path(r2), path(host_index), path(hbv_index)

    output:
    tuple val(meta), path("${meta.id}_host.sam"), path("${meta.id}_hbv.sam"), emit: aligned_sam
    tuple val(meta), path("${meta.id}_host.bam"), path("${meta.id}_hbv.bam"), emit: aligned_bam
    tuple val(meta), path("${meta.id}_first_align.log"), emit: first_align_log

    script:
    """
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${host_index}/${params.host_bwa_prefix} ${r1} ${r2} 2> ${meta.id}_host.log > ${meta.id}_host.sam
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${hbv_index}/genome.fa ${r1} ${r2} 2> ${meta.id}_hbv.log > ${meta.id}_hbv.sam
    samtools view -bS ${meta.id}_host.sam > ${meta.id}_host.bam
    samtools view -bS ${meta.id}_hbv.sam > ${meta.id}_hbv.bam
    cat ${meta.id}_host.log ${meta.id}_hbv.log > ${meta.id}_first_align.log
    """
}


process FILTER_CHIMERIC {
    tag "Filtering chimeric reads of ${meta.id}"

    input:
    tuple val(meta), path(host_sam), path(hbv_sam), path(trim_r1), path(trim_r2)

    output:
    tuple val(meta), path("${meta.id}_chimeric_R1.fq.gz"), path("${meta.id}_chimeric_R2.fq.gz"), emit: chimeric_fq
    tuple val(meta), path("${meta.id}_chimeric.log"), emit: chimeric_log
    tuple val(meta), path("${meta.id}_molecule_map.tsv"), emit: molecule_map

    script:
    """
    filter_chimeric.py ${host_sam} ${hbv_sam} --r1 ${trim_r1} --r2 ${trim_r2} 2> ${meta.id}_chimeric.log
    mv molecule_map.tsv ${meta.id}_molecule_map.tsv
    mv chimeric_R1.fq.gz ${meta.id}_chimeric_R1.fq.gz
    mv chimeric_R2.fq.gz ${meta.id}_chimeric_R2.fq.gz
    """
}


process QC_CAPTURE {
    tag "HBV-capture QC of ${meta.id}"

    input:
    tuple val(meta), path(trim_R1), path(hbv_sam)

    output:
    tuple val(meta), path("${meta.id}_hbv_qc.tsv"), emit: qc_tsv

    script:
    """
    qc_capture.py ${meta.id} ${trim_R1} ${hbv_sam} > ${meta.id}_hbv_qc.tsv
    """
}


process PE_ASSEMBLE {
    tag "PE-assembling chimeric reads of ${meta.id} with bbmerge"

    input:
    tuple val(meta), path(chim_R1), path(chim_R2)

    output:
    tuple val(meta), path("${meta.id}_merged.fq.gz"), emit: merged_fq
    tuple val(meta), path("${meta.id}_unmerged_R1.fq.gz"), path("${meta.id}_unmerged_R2.fq.gz"), emit: unmerged_fq

    tuple val(meta), path("${meta.id}_assembly.log"), path("${meta.id}_assembly_qc.json"), emit: assembly_qc

    script:
    """
    assemble_pairs.py ${chim_R1} ${chim_R2} \\
        --prefix ${meta.id} \\
        --min-overlap ${params.minOverlap} \\
        --max-mismatch ${params.maxMismatchRate} \\
        --threads ${task.cpus} \\
        --keep-unmerged ${params.keepUnmerged}
    """
}


process PE_ASSEMBLE_CUSTOM {
    tag "PE-assembling chimeric reads of ${meta.id} (custom method)"

    input:
    tuple val(meta), path(chim_R1), path(chim_R2)

    output:
    tuple val(meta), path("${meta.id}_merged.fq.gz"), emit: merged_fq
    tuple val(meta), path("${meta.id}_unmerged_R1.fq.gz"), path("${meta.id}_unmerged_R2.fq.gz"), emit: unmerged_fq

    tuple val(meta), path("${meta.id}_assembly.log"), path("${meta.id}_assembly_qc.json"), emit: assembly_qc

    script:
    """
    pe_assemble.py ${chim_R1} ${chim_R2} \\
        --prefix ${meta.id} \\
        --min-overlap ${params.minOverlap} \\
        --max-mismatch ${params.maxMismatchRate} \\
        --keep-unmerged ${params.keepUnmerged}
    """
}


process REMAP {
    tag "Re-mapping merged reads of ${meta.id} and calling breakpoints"

    input:
    tuple val(meta), path(merged_fq), path(host_index), path(hbv_index)

    output:
    tuple val(meta), path("${meta.id}_remap_host.sam"), path("${meta.id}_remap_hbv.sam"), emit: remap_sam
    tuple val(meta), path("${meta.id}_breakpoints.txt"), emit: breakpoints

    script:
    """
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${host_index}/${params.host_bwa_prefix} ${merged_fq} 2> ${meta.id}_remap_host.log > ${meta.id}_remap_host.sam
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${hbv_index}/genome.fa ${merged_fq} 2> ${meta.id}_remap_hbv.log > ${meta.id}_remap_hbv.sam
    call_breakpoints.py ${meta.id}_remap_host.sam ${meta.id}_remap_hbv.sam --min-match ${params.minMatch} --host-index ${host_index}/${params.host_bwa_prefix} --min-score ${params.bwaMinScore} --sample ${meta.id} --min-unique ${params.minUnique} --max-overlap ${params.maxJunctionOverlap} --score-margin ${params.junctionScoreMargin} --min-mapq ${params.minMapq} --alt-score-delta ${params.altScoreDelta} > ${meta.id}_breakpoints.txt
    """
}


process MERGE_SIGNAL {
    tag "Merging breakpoints within ${params.mergeDist} bp of ${meta.id}"

    input:
    tuple val(meta), path(bp)

    output:
    tuple val(meta), path("${meta.id}_breakpoints_merged.txt"), emit: merged_breakpoints

    script:
    """
    cat ${bp} | merge_breakpoints.py --dist ${params.mergeDist} > ${meta.id}_breakpoints_merged.txt
    """
}


process COUNT_UNMERGED_NSS {
    tag "Merging unmerged chimeric reads into breakpoints/NSS of ${meta.id}"

    input:
    tuple val(meta), path(unmerged_R1), path(unmerged_R2), path(host_index), path(hbv_index), path(breakpoints)

    output:
    tuple val(meta), path("${meta.id}_support_host.sam"), path("${meta.id}_support_hbv.sam"), emit: support_sam
    tuple val(meta), path("${meta.id}_breakpoints_combined.txt"), emit: breakpoints

    script:
    """
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${host_index}/${params.host_bwa_prefix} ${unmerged_R1} ${unmerged_R2} 2> ${meta.id}_support_host.log > ${meta.id}_support_host.sam
    bwa mem -SP -Y -T ${params.bwaMinScore} -t ${task.cpus} ${hbv_index}/genome.fa ${unmerged_R1} ${unmerged_R2} 2> ${meta.id}_support_hbv.log > ${meta.id}_support_hbv.sam
    count_support.py ${meta.id}_support_host.sam ${meta.id}_support_hbv.sam ${breakpoints} --min-match ${params.minMatch} --dist ${params.supportDist} --merge-dist ${params.mergeDist} --host-index ${host_index}/${params.host_bwa_prefix} --min-score ${params.bwaMinScore} --sample ${meta.id} --min-unique ${params.minUnique} --max-overlap ${params.maxJunctionOverlap} --score-margin ${params.junctionScoreMargin} --min-mapq ${params.minMapq} --alt-score-delta ${params.altScoreDelta} > ${meta.id}_breakpoints_combined.txt
    """
}


process FILTER_SIGNAL {
    tag "Normalizing NSS to NNSS of ${meta.id}"

    input:
    tuple val(meta), path(bp_combined), path(trim_R1)

    output:
    tuple val(meta), path("${meta.id}_breakpoints_final.txt"), emit: final_breakpoints
    tuple val(meta), path("${meta.id}_breakpoints_primary.txt"), emit: primary_breakpoints
    tuple val(meta), path("${meta.id}_nnss_summary.txt"), emit: nnss_summary

    script:
    """
    EFR=\$(zcat ${trim_R1} | wc -l)
    EFR=\$(( EFR / 4 ))
    cat ${bp_combined} | filter_signal.py \${EFR} > ${meta.id}_breakpoints_final.txt
    collapse_secondary.py ${meta.id}_breakpoints_final.txt > ${meta.id}_breakpoints_primary.txt
    printf '%s\\t%s\\n' ${meta.id} \${EFR} > ${meta.id}_nnss_summary.txt
    """
}




workflow {

    main:
    log.info """\
      nftide-hivid (HBV integration breakpoint detection)
      ===================================
      projectDir             :  ${projectDir}
      workingDir             :  ${workflow.outputDir}
      input csv              :  ${params.input_csv}
      HBV meta csv           :  ${params.HBVfa_csv}
      host bwa index         :  ${params.host_bwa_dir}/${params.host_bwa_prefix}
      bwa min score (-T)     :  ${params.bwaMinScore}
      min match (bp)         :  ${params.minMatch}
      PE min overlap         :  ${params.minOverlap}
      PE max mismatch rate   :  ${params.maxMismatchRate}
      merge distance (bp)    :  ${params.mergeDist}
      NNSS threshold         :  ${params.nnssThreshold}  (deprecated - final table is not filtered)
      PE assembly method     :  ${params.peAssembleCustom ? 'custom script' : 'bbmerge.sh'}
      keep unmerged reads    :  ${params.keepUnmerged}
    """.stripIndent()

    ch_read_pairs = channel.fromPath(params.input_csv)
    .splitCsv(header:true)
    .map { row -> 
        [
            row.sample,
            row
        ]
    }
    .groupTuple()
    .map { _sample, rows -> 
        rows.withIndex().collect { row, index ->
            row + [rep: index + 1]
        }
    }
    .flatMap { item -> item }
   .map { row -> 

        [
            [
                id: row.sample,
                rep: row.rep,

            ], 
            [
                file(row.fastq_1, checkIfExists: true), 
                file(row.fastq_2, checkIfExists: true)
            ]
        ]
    }
    .map{meta, files -> [meta.subMap(['id']), files]}
    .groupTuple()
    .map { meta, filePairs ->
        [ meta, filePairs.collect { pair -> pair[0] }, filePairs.collect { pair -> pair[1] }]
    }

    MERGE_FQ(ch_read_pairs)
    CUTADAPT(MERGE_FQ.out.merged_fq)

    // HBV reference genome per sample (from meta_all_assembled_fa.csv)
    ch_hbv = channel.fromPath(params.HBVfa_csv)
        .splitCsv(header: true)
        .map { row -> [ [id: row.sample], row.sequence ] }

    // sanity check: every sample must have an HBV reference genome
    // (read both CSVs directly; simple comma-separated format)
    def _samples = new File(params.input_csv).readLines().drop(1).collect { line -> line.split(',')[0].trim() }
    def _hbv_samples = new File(params.HBVfa_csv).readLines().drop(1).collect { line -> line.split(',')[0].trim() }
    def missing_ids = _samples - _hbv_samples
    CHECK_SAMPLES(missing_ids ? missing_ids.join(' ') : '')

    // step 0 (helper): extract per-sample HBV genome and build HBV index
    // (the host bwa index is prebuilt and shared by all samples)
    EXTRACT_HBV_FA(ch_hbv)
    BUILD_HBV_INDEX(EXTRACT_HBV_FA.out.hbv_fa)

    ch_host_idx = channel.value(file(params.host_bwa_dir))
    ch_hbv_idx = BUILD_HBV_INDEX.out.hbv_index

    ch_trimmed = CUTADAPT.out.trimmed_reads

    // step 2a: first-time alignment (bwa mem on host + HBV)
    // note: cross() silently emits nothing in some Nextflow versions -> use combine()
    ch_first_in = ch_trimmed.join(ch_hbv_idx)
        .combine(ch_host_idx)
        .map { it ->
            def (meta, r1, r2, hbv_idx, host_idx) = it
            [ meta, r1, r2, host_idx, hbv_idx ]
        }
    FIRST_ALIGN(ch_first_in)

    // step 2b: filter chimeric read pairs (host + HBV)
    // (separate process so changes here do NOT re-run the expensive bwa step)
    FILTER_CHIMERIC(FIRST_ALIGN.out.aligned_sam.join(ch_trimmed))

    // step 2c: HBV-capture QC (total read pairs, HBV-aligned pairs, proportion)
    ch_qc_in = FIRST_ALIGN.out.aligned_sam
        .map { meta, _host_sam, hbv_sam -> [ meta, hbv_sam ] }
        .join(ch_trimmed.map { meta, r1, _r2 -> [ meta, r1 ] })
        .map { it ->
            def (meta, hbv_sam, r1) = it
            [ meta, r1, hbv_sam ]
        }
    QC_CAPTURE(ch_qc_in)

    // step 3: paired-end assembly of chimeric reads
    // default: bbmerge.sh; only use the custom script when peAssembleCustom = true
    // Sequenced junctions from merged reads or individual mates take priority.
    // Discordant pairs support compatible sequenced junctions or retain estimated coordinates.
    if (params.peAssembleCustom) {
        PE_ASSEMBLE_CUSTOM(FILTER_CHIMERIC.out.chimeric_fq)
        ch_merged = PE_ASSEMBLE_CUSTOM.out.merged_fq
        ch_unmerged = PE_ASSEMBLE_CUSTOM.out.unmerged_fq
        ch_assembly_qc = PE_ASSEMBLE_CUSTOM.out.assembly_qc
    } else {
        PE_ASSEMBLE(FILTER_CHIMERIC.out.chimeric_fq)
        ch_merged = PE_ASSEMBLE.out.merged_fq
        ch_unmerged = PE_ASSEMBLE.out.unmerged_fq
        ch_assembly_qc = PE_ASSEMBLE.out.assembly_qc
    }

    // step 4: re-mapping + breakpoint calling (on merged reads only)
    ch_remap_in = ch_merged.join(ch_hbv_idx)
        .combine(ch_host_idx)
        .map { it ->
            def (meta, merged_fq, hbv_idx, host_idx) = it
            [ meta, merged_fq, host_idx, hbv_idx ]
        }
    REMAP(ch_remap_in)

    // step 5b: count unmerged chimeric pairs as NSS support for breakpoints.
    // Uses the RAW breakpoints from REMAP: all breakpoints (merged-read and
    // unmerged-defined) are called first, then merged together in step 5.
    ch_support_in = ch_unmerged.join(ch_hbv_idx)
        .combine(ch_host_idx)
        .map { it ->
            def (meta, u1, u2, hbv_idx, host_idx) = it
            [ meta, u1, u2, host_idx, hbv_idx ]
        }
        .join(REMAP.out.breakpoints.map { meta, bp -> [ meta, bp ] })
    COUNT_UNMERGED_NSS(ch_support_in)

    // step 5: merge ALL breakpoints (merged + unmerged) within 20 bp;
    // both references and orientations must agree; NSS counts molecule proxies
    MERGE_SIGNAL(COUNT_UNMERGED_NSS.out.breakpoints)

    // step 6: NNSS normalization (NSS already includes unmerged support)
    ch_efr = ch_trimmed.map { meta, r1, _r2 -> [ meta, r1 ] }
    ch_filter_in = MERGE_SIGNAL.out.merged_breakpoints.join(ch_efr)
    FILTER_SIGNAL(ch_filter_in)



    publish:
    out_trimmed_fastqs = CUTADAPT.out.trimmed_reads
    out_cutadapt_logs = CUTADAPT.out.cutadapt_log
    out_first_align_logs = FIRST_ALIGN.out.first_align_log
    out_first_align_bam = FIRST_ALIGN.out.aligned_bam
    out_hbv_qc = QC_CAPTURE.out.qc_tsv
    out_chimeric_fastqs = FILTER_CHIMERIC.out.chimeric_fq
    out_chimeric_logs = FILTER_CHIMERIC.out.chimeric_log
    out_molecule_map = FILTER_CHIMERIC.out.molecule_map
    out_assembly_qc = ch_assembly_qc
    out_assembled_reads = ch_merged
    out_unmerged_reads = ch_unmerged
    out_hbv_fa = EXTRACT_HBV_FA.out.hbv_fa
    out_remap_sam = REMAP.out.remap_sam
    out_support_sam = COUNT_UNMERGED_NSS.out.support_sam
    out_raw_breakpoints = REMAP.out.breakpoints
    out_merged_breakpoints = MERGE_SIGNAL.out.merged_breakpoints
    out_combined_breakpoints = COUNT_UNMERGED_NSS.out.breakpoints
    out_final_breakpoints = FILTER_SIGNAL.out.final_breakpoints
    out_primary_breakpoints = FILTER_SIGNAL.out.primary_breakpoints
    out_nnss_summary = FILTER_SIGNAL.out.nnss_summary


}

output {
    out_assembly_qc {
        path { meta, _log, _qc -> "${meta.id}/assembly" }
    }
    out_molecule_map {
        path { meta, _f -> "${meta.id}/qc" }
    }
    out_trimmed_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_cutadapt_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_first_align_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_first_align_bam {
        path { meta, _f1, _f2 -> "${meta.id}/sam" }
    }
    out_hbv_qc {
        path { meta, _f -> "${meta.id}/qc" }
    }
    out_chimeric_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_chimeric_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_assembled_reads {
        path { meta, _f -> "${meta.id}/assembly" }
    }
    out_unmerged_reads {
        path { meta, _f1, _f2 -> "${meta.id}/assembly" }
    }
    out_hbv_fa {
        path { meta, _f -> _f >> "${meta.id}/fastqs/${meta.id}_hbv.fa" }
    }
    out_remap_sam {
        path { meta, _f1, _f2 -> "${meta.id}/sam" }
    }
    out_support_sam {
        path { meta, _f1, _f2 -> "${meta.id}/sam" }
    }
    out_raw_breakpoints {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
    out_merged_breakpoints {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
    out_combined_breakpoints {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
    out_final_breakpoints {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
    out_primary_breakpoints {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
    out_nnss_summary {
        path { meta, _f -> "${meta.id}/breakpoints" }
    }
}