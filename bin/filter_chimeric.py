#!/usr/bin/env python3
"""Extract chimeric read pairs (host + HBV) from two BWA-MEM SAM files.

Reads paired-end SAM alignments against the host genome (host.sam) and the
HBV genome (hbv.sam) and writes out, as gzipped fastq, every read pair that
aligns to BOTH genomes (chimeric = viral integration junction candidates).

This corresponds to the HIVID "raw-mapping" step: keep chimeric paired-end
reads whose sequence aligns partially to the human genome and partially to
the HBV genome.
"""
import gzip
import sys


def collect(path):
    """Return {qname: {'r1': [rec], 'r2': [rec]}} keeping the first record per segment."""
    recs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            qname = f[0]
            flag = int(f[1])
            rname = f[2]
            if rname == '*':
                continue
            if flag & 0x40:
                seg = 'r1'
            elif flag & 0x80:
                seg = 'r2'
            else:
                seg = 's'
            recs.setdefault(qname, {}).setdefault(seg, []).append(f)
    return recs


def main():
    if len(sys.argv) < 3:
        sys.stderr.write('usage: filter_chimeric.py host.sam hbv.sam\n')
        sys.exit(1)
    host = collect(sys.argv[1])
    hbv = collect(sys.argv[2])

    n_chimeric = 0
    with gzip.open('chimeric_R1.fq.gz', 'wt') as out1, \
         gzip.open('chimeric_R2.fq.gz', 'wt') as out2:
        for qname in host:
            if qname not in hbv:
                continue
            n_chimeric += 1
            r1 = host[qname].get('r1') or hbv[qname].get('r1')
            r2 = host[qname].get('r2') or hbv[qname].get('r2')
            if r1:
                rec = r1[0]
                out1.write('@{0}\n{1}\n+\n{2}\n'.format(qname, rec[9], rec[10]))
            if r2:
                rec = r2[0]
                out2.write('@{0}\n{1}\n+\n{2}\n'.format(qname, rec[9], rec[10]))
    sys.stderr.write('chimeric pairs: {0}\n'.format(n_chimeric))


if __name__ == '__main__':
    main()
