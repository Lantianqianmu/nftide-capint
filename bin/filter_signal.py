#!/usr/bin/env python3
"""Normalize NSS to NNSS (HIVID step 6).

    NNSS = NSS * 10^6 / EFR

EFR is the number of effective read pairs (read pairs used to search for
breakpoints). This script reports ALL input breakpoints WITHOUT filtering;
it only appends the EFR and NNSS columns.

Input columns (from count_support.py, combined table):
    host_chr  host_pos  host_strand  hbv_pos  hbv_strand  nss
    source  host_mapq  hbv_mapq  event_id  event_type
Output adds EFR and NNSS columns at the end.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('efr', type=int, help='effective read pairs')
    p.add_argument('--nnss-threshold', type=float, default=1.0,
                   help='(deprecated) NNSS threshold; no longer filters, '
                        'kept for backward compatibility')
    args = p.parse_args()

    out = sys.stdout
    out.write('host_chr\thost_pos\thost_strand\thbv_pos\thbv_strand'
              '\tnss\tsource\thost_mapq\thbv_mapq'
              '\tevent_id\tevent_type\tefr\tnnss\n')
    for line in sys.stdin:
        line = line.rstrip('\n')
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 7:
            continue
        if f[0] == 'host_chr':
            continue  # skip input header
        host_chr, host_pos, hstrand, hbv_pos, vstrand, nss = f[:6]
        source = f[6] if len(f) > 6 else 'merged'
        host_mapq = int(f[7]) if len(f) > 7 and f[7] else 0
        hbv_mapq = int(f[8]) if len(f) > 8 and f[8] else 0
        event_id = f[9] if len(f) > 9 else ''
        event_type = f[10] if len(f) > 10 else 'primary'
        try:
            nss_v = int(nss)
        except ValueError:
            continue
        nnss = nss_v * 1e6 / float(args.efr) if args.efr else 0.0
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12:.3f}\n'.format(
            host_chr, host_pos, hstrand, hbv_pos, vstrand,
            nss_v, source, host_mapq, hbv_mapq, event_id, event_type,
            args.efr, nnss))


if __name__ == '__main__':
    main()
