#!/usr/bin/env python3
"""Report trimmed input, HBV-aligned and sequence-deduplicated HBV pair counts.

Unique means exact full-pair sequence proxies, not unique mapping or UMI counts.
The first four output columns retain their existing order.
"""
import argparse
import sys
from assemble_pairs import paired_records
from filter_chimeric import fastq, mapped_names
from junctions import molecule_id


def capture_metrics(r1, r2, sam):
    hbv_names = mapped_names(sam)
    total = 0
    molecules = set()
    for first, second in paired_records(fastq(r1), fastq(r2)):
        total += 1
        if first[0] in hbv_names:
            molecules.add(molecule_id('', [first[1], second[1]]))
    return total, len(hbv_names), len(molecules)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('sample')
    p.add_argument('trim_r1')
    p.add_argument('hbv_sam')
    p.add_argument('--trim-r2', required=True,
                   help='matching trimmed R2 FASTQ; required for full-pair deduplication')
    args = p.parse_args()
    total, hbv, unique = capture_metrics(args.trim_r1, args.trim_r2, args.hbv_sam)
    prop = hbv / total if total else 0.0
    sys.stdout.write('sample\ttotal_read_pairs\thbv_read_pairs\tprop_hbv_reads\tunique_hbv_read_pairs\n')
    sys.stdout.write(f'{args.sample}\t{total}\t{hbv}\t{prop:.6f}\t{unique}\n')


if __name__ == '__main__':
    main()
