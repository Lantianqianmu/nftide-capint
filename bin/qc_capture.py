#!/usr/bin/env python3
"""Capture basic HBV-capture QC metrics after the first alignment.

Reports, as a single-row TSV with a header:
    sample            sample id
    total_read_pairs  number of read pairs in the trimmed R1 fastq
    hbv_read_pairs    number of read pairs with at least one read aligned to HBV
    prop_hbv_reads    hbv_read_pairs / total_read_pairs
"""
import argparse
import gzip
import sys


def count_read_pairs(fastq):
    """Number of read pairs = number of records in R1 fastq."""
    n = 0
    opener = gzip.open if fastq.endswith('.gz') else open
    with opener(fastq, 'rt') as fh:
        for _line in fh:
            n += 1
    return n // 4


def count_hbv_pairs(sam):
    """Read pairs (unique qnames) with at least one mapped HBV alignment."""
    seen = set()
    with open(sam) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6 or f[2] == '*' or f[5] == '*':
                continue  # unmapped
            seen.add(f[0])
    return len(seen)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('sample')
    p.add_argument('trim_r1')
    p.add_argument('hbv_sam')
    args = p.parse_args()

    total = count_read_pairs(args.trim_r1)
    hbv = count_hbv_pairs(args.hbv_sam)
    prop = hbv / total if total else 0.0

    sys.stdout.write('sample\ttotal_read_pairs\thbv_read_pairs\tprop_hbv_reads\n')
    sys.stdout.write('{0}\t{1}\t{2}\t{3:.6f}\n'.format(
        args.sample, total, hbv, prop))


if __name__ == '__main__':
    main()
