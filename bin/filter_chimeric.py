#!/usr/bin/env python3
"""Select mapped-to-both candidate pairs; recover original FASTQs and deduplicate.

Exact full-pair sequence duplicates are collapsed before assembly. This is a
conservative sequence-based molecule proxy; no UMI information is inferred.
"""
import argparse
import gzip
import itertools
import sys
from junctions import molecule_id, rc


def mapped_names(path):
    names = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            if len(f) >= 11 and not int(f[1]) & 4 and f[2] != '*' and f[5] != '*':
                names.add(f[0])
    return names


def fastq(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        while True:
            name = handle.readline().rstrip()
            if not name:
                return
            seq = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            qual = handle.readline().rstrip()
            if not name.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Invalid FASTQ record in ' + str(path))
            name = name[1:].split()[0]
            if name.endswith(('/1','/2')):
                name = name[:-2]
            yield name, seq, qual


def recovered(paths, candidates):
    # Backward-compatible SAM-only use. Ignore hard-clipped / alternative
    # records so only complete primary sequences are used for reconstruction.
    reads = {}
    for path in paths:
        with open(path) as handle:
            for line in handle:
                if line.startswith('@'):
                    continue
                f = line.rstrip().split('\t')
                if len(f) < 11 or f[0] not in candidates:
                    continue
                flag = int(f[1])
                if flag & (256|2048) or 'H' in f[5] or f[9] == '*' or f[10] == '*':
                    continue
                if not flag & (64|128):
                    continue
                seq, qual = f[9], f[10]
                if flag & 16:
                    seq, qual = rc(seq), qual[::-1]
                reads.setdefault(f[0], {})[1 if flag & 128 else 0] = (f[0], seq, qual)
    for name in sorted(reads):
        if len(reads[name]) == 2:
            yield reads[name][0], reads[name][1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host_sam'); p.add_argument('hbv_sam')
    p.add_argument('--r1'); p.add_argument('--r2')
    args = p.parse_args()
    if bool(args.r1) != bool(args.r2):
        p.error('--r1 and --r2 must be supplied together')
    candidates = mapped_names(args.host_sam) & mapped_names(args.hbv_sam)
    pairs = itertools.zip_longest(fastq(args.r1), fastq(args.r2)) if args.r1 else recovered((args.host_sam,args.hbv_sam), candidates)
    unique = {}
    with open('molecule_map.tsv', 'w') as mapping:
        mapping.write('read_id\tmolecule_id\n')
        for a, b in pairs:
            if a is None or b is None or a[0] != b[0]:
                raise ValueError('FASTQ mates are not synchronized')
            if a[0] not in candidates:
                continue
            mid = molecule_id('', [a[1],b[1]])
            mapping.write(a[0]+'\t'+mid+'\n')
            if mid in unique:
                unique[mid][2] += 1
            else:
                unique[mid] = [a,b,1]
    with gzip.open('chimeric_R1.fq.gz','wt') as out1, gzip.open('chimeric_R2.fq.gz','wt') as out2:
        for mid, (a,b,count) in unique.items():
            name = a[0]+'|mol='+mid+'|copies='+str(count)
            for rec, out in ((a,out1),(b,out2)):
                out.write('@'+name+'\n'+rec[1]+'\n+\n'+rec[2]+'\n')
    sys.stderr.write('candidate pairs: %d; sequence-deduplicated molecules: %d\n' % (sum(x[2] for x in unique.values()),len(unique)))

if __name__ == '__main__':
    main()
