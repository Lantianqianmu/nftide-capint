#!/usr/bin/env python3
"""PE-assembly of chimeric read pairs following the HIVID procedure.

Both ends are put onto the same strand (R2 is reverse-complemented). If the
tail of the upstream end overlaps the head of the downstream end by more than
5 bp with a mismatch rate < 0.2 the two ends are spliced into one continuous
sequence, the PE-assembled read. Pairs that cannot be spliced are written out
as two separate reads so that the breakpoint search is not lost.

Output: a single gzipped fastq where
  - spliced reads are written as one record        (name)
  - non-spliced pairs are written as two records   (name/1 and name/2)
"""
import argparse
import gzip
import sys

COMP = str.maketrans('ACGTNacgtnRYSWKMryswkmBDHVbdhv',
                     'TGCANtgcanYRSWMKyrswmkVHDBvhdb')


def revcomp(s):
    return s.translate(COMP)[::-1]


def merge(a, b, min_overlap, max_mismatch):
    """Splice b onto the tail of a if the tail of a overlaps the head of b."""
    max_ov = min(len(a), len(b))
    for ov in range(max_ov, min_overlap - 1, -1):
        mism = sum(1 for x, y in zip(a[-ov:], b[:ov]) if x != y)
        if mism / float(ov) <= max_mismatch:
            return a + b[ov:]
    return None


def assemble(r1, r2, min_overlap, max_mismatch):
    r2rc = revcomp(r2)
    m = merge(r1, r2rc, min_overlap, max_mismatch)
    if m:
        return [m]
    r1rc = revcomp(r1)
    m2 = merge(r2, r1rc, min_overlap, max_mismatch)
    if m2:
        return [m2]
    return [r1, r2]


def read_fastq(fh):
    while True:
        name = fh.readline()
        if not name:
            return
        seq = fh.readline().strip()
        fh.readline()          # '+'
        qual = fh.readline().strip()
        yield name.strip(), seq, qual


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('r1', help='chimeric R1 fastq.gz')
    p.add_argument('r2', help='chimeric R2 fastq.gz')
    p.add_argument('--out', default='merged.fq.gz',
                   help='output of spliced (PE-assembled) reads')
    p.add_argument('--out-u1', default='unmerged_R1.fq.gz',
                   help='output of non-overlapping pair R1 (when --keep-unmerged)')
    p.add_argument('--out-u2', default='unmerged_R2.fq.gz',
                   help='output of non-overlapping pair R2 (when --keep-unmerged)')
    p.add_argument('--min-overlap', type=int, default=6,
                   help='minimum overlap (>5 bp) required for splicing')
    p.add_argument('--max-mismatch', type=float, default=0.2,
                   help='maximum mismatch rate allowed in the overlap')
    p.add_argument('--keep-unmerged', action='store_true',
                   help='also write non-overlapping pairs as /1 /2 records '
                        '(default: only spliced reads are written)')
    args = p.parse_args()

    recs = {}
    with gzip.open(args.r1, 'rt') as fh:
        for name, seq, qual in read_fastq(fh):
            base = name[:-2] if name.endswith('/1') else name
            recs.setdefault(base, {})['r1'] = seq
    with gzip.open(args.r2, 'rt') as fh:
        for name, seq, qual in read_fastq(fh):
            base = name[:-2] if name.endswith('/2') else name
            recs.setdefault(base, {})['r2'] = seq

    n_pairs = 0
    n_spliced = 0
    out_u1 = gzip.open(args.out_u1, 'wt')   # always created (empty if not kept)
    out_u2 = gzip.open(args.out_u2, 'wt')
    with gzip.open(args.out, 'wt') as out:
        for base, pair in recs.items():
            if 'r1' not in pair or 'r2' not in pair:
                continue
            n_pairs += 1
            parts = assemble(pair['r1'], pair['r2'],
                             args.min_overlap, args.max_mismatch)
            if len(parts) == 1:
                n_spliced += 1
                out.write('@{0}\n{1}\n+\n{2}\n'.format(
                    base, parts[0], 'I' * len(parts[0])))
            elif args.keep_unmerged:
                out_u1.write('@{0}/1\n{1}\n+\n{2}\n'.format(
                    base, parts[0], 'I' * len(parts[0])))
                out_u2.write('@{0}/2\n{1}\n+\n{2}\n'.format(
                    base, parts[1], 'I' * len(parts[1])))
    out_u1.close()
    out_u2.close()
    sys.stderr.write('chimeric pairs: {0}, spliced: {1}\n'.format(n_pairs, n_spliced))


if __name__ == '__main__':
    main()
