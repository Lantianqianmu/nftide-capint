#!/usr/bin/env python3
"""Shared helper: re-map reads to the host with `bwa mem -a` to enumerate ALL
possible alignment positions.

Used for secondary-event calling: when a read's primary host alignment has
MAPQ == 0 (ambiguous / multi-mapping, AS == XS), we re-map the read and report
every other locus it hits as a "secondary" breakpoint candidate.
"""
import os
import re
import subprocess
import tempfile

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def cigar_matches(cigar):
    """Number of reference-matching bases (M, =, X) of a CIGAR string."""
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in 'M=X')


def remap_host_all(reads, host_index, min_score, threads=2):
    """Map every read to the host with `bwa mem -a` and return ALL alignments.

    reads      : list of (read_id, sequence)
    host_index : path to the host bwa index prefix (e.g. /path/hg38.fa)
    min_score  : bwa mem -T score threshold (only report alignments >= score)

    Returns {read_id: [(is_secondary, rname, pos, strand, mapq, ref_len)]}
    where is_secondary=True for alignments flagged 0x100 (suboptimal hits;
    the primary alignment is the one without that flag). Supplementary
    alignments (flag 2048) are ignored.
    """
    if not reads or not host_index:
        return {}
    fa = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False)
    try:
        for rid, seq in reads:
            fa.write('>{0}\n{1}\n'.format(rid, seq))
        fa.flush()
        cmd = ['bwa', 'mem', '-a', '-T', str(min_score), '-t', str(threads),
               host_index, fa.name]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
    finally:
        fa.close()
        if os.path.exists(fa.name):
            os.unlink(fa.name)

    out = {}
    for line in proc.stdout.decode('utf-8', 'replace').splitlines():
        if line.startswith('@'):
            continue
        f = line.split('\t')
        if len(f) < 6:
            continue
        flag = int(f[1])
        if flag & 0x800:          # supplementary alignment -> ignore
            continue
        rname = f[2]
        cigar = f[5]
        if rname == '*' or cigar == '*':
            continue
        ref_len = cigar_matches(cigar)
        if ref_len <= 0:
            continue
        out.setdefault(f[0], []).append(
            (bool(flag & 0x100), rname, int(f[3]),
             '-' if flag & 0x10 else '+', int(f[4]), ref_len))
    return out
