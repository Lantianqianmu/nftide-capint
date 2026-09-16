#!/usr/bin/env python3
"""Custom ungapped suffix/prefix paired-read assembler.

R2 is reverse-complemented; the longest overlap meeting min-overlap and
max-mismatch is used. This simple algorithm is not BBMerge. Output names,
molecule tags, mate suffixes and QC follow the Nextflow assembly contract.
"""
import argparse
import gzip
import json
from pathlib import Path
from assemble_pairs import paired_records, write_record
from filter_chimeric import fastq

COMP = str.maketrans('ACGTNacgtnRYSWKMryswkmBDHVbdhv',
                     'TGCANtgcanYRSWMKyrswmkVHDBvhdb')


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def overlap(a, b, minimum, maximum):
    for size in range(min(len(a), len(b)), minimum - 1, -1):
        # Ambiguous bases do not provide evidence of an overlap.
        mismatches = sum(x != y or x not in 'ACGT'
                         for x, y in zip(a[-size:].upper(), b[:size].upper()))
        if mismatches / float(size) <= maximum:
            return size
    return None


def consensus(first, second, minimum, maximum):
    name, a, aq = first
    _, raw_b, raw_bq = second
    b, bq = revcomp(raw_b), raw_bq[::-1]
    size = overlap(a, b, minimum, maximum)
    if size is None:
        return None
    seq, qual = list(a[:-size]), list(aq[:-size])
    for x, y, qx, qy in zip(a[-size:], b[:size], aq[-size:], bq[:size]):
        if x.upper() == y.upper():
            seq.append(x)
            qual.append(max(qx, qy))
        else:
            # Prefer the higher-quality observation. Equal-quality conflicts
            # remain uncertain, rather than manufacturing a high-quality base.
            seq.append(x if qx > qy else y if qy > qx else 'N')
            qual.append(chr(33 + abs(ord(qx) - ord(qy))))
    return name, ''.join(seq) + b[size:], ''.join(qual) + bq[size:]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('r1'); p.add_argument('r2')
    p.add_argument('--prefix', help='pipeline output prefix, including log and QC')
    p.add_argument('--out', default='merged.fq.gz')
    p.add_argument('--out-u1', default='unmerged_R1.fq.gz')
    p.add_argument('--out-u2', default='unmerged_R2.fq.gz')
    p.add_argument('--min-overlap', type=int, default=6)
    p.add_argument('--max-mismatch', type=float, default=.2)
    p.add_argument('--keep-unmerged', nargs='?', const='true',
                   choices=('true', 'false'), default='false')
    args = p.parse_args()
    if args.min_overlap < 1 or not 0 <= args.max_mismatch <= 1:
        p.error('min-overlap must be positive; max-mismatch must be between 0 and 1')
    outputs = ([args.prefix + suffix for suffix in
                ('_merged.fq.gz', '_unmerged_R1.fq.gz', '_unmerged_R2.fq.gz')]
               if args.prefix else [args.out, args.out_u1, args.out_u2])
    temporary = [Path(path + '.tmp') for path in outputs]
    keep = args.keep_unmerged == 'true'
    total = merged = 0
    try:
        with gzip.open(temporary[0], 'wt') as out, \
             gzip.open(temporary[1], 'wt') as u1, \
             gzip.open(temporary[2], 'wt') as u2:
            for first, second in paired_records(fastq(args.r1), fastq(args.r2)):
                total += 1
                result = consensus(first, second, args.min_overlap, args.max_mismatch)
                if result is not None:
                    merged += 1
                    write_record(out, result)
                elif keep:
                    write_record(u1, first, '/1')
                    write_record(u2, second, '/2')
        for source, destination in zip(temporary, outputs):
            source.replace(destination)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
    summary = dict(input_pairs=total, merged_pairs=merged,
                   unmerged_pairs=total-merged, keep_unmerged=keep,
                   emitted_unmerged_pairs=total-merged if keep else 0)
    report = json.dumps(summary, indent=2) + '\n'
    if args.prefix:
        Path(args.prefix + '_assembly_qc.json').write_text(report)
        Path(args.prefix + '_assembly.log').write_text('Custom ungapped overlap assembler\n' + report)
    print(report, end='')


if __name__ == '__main__':
    main()
