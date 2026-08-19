#!/usr/bin/env python3
"""Call HBV integration breakpoints from BWA-MEM alignments (HIVID step 4).

Reads SAM alignments of the assembled reads against the host (host.sam) and
the HBV (hbv.sam) genomes. A breakpoint is reported for every assembled read
whose alignment covers at least --min-match bp on BOTH genomes. The joint
position of the human and the HBV sequence is the reported breakpoint.

Supplementary alignments (SAM flag 2048) are ignored: only the primary
alignment of each read is considered.

SECONDARY EVENTS: when a called breakpoint has host_mapq == 0 (ambiguous /
multi-mapping host alignment, AS == XS) the supporting read is re-mapped to
the host with `bwa mem -a` and every additional locus it hits (>= --min-match
bp, flagged 0x100) is reported as a separate "secondary" breakpoint row with
event_type = secondary. This requires --host-index; without it only primary
rows are emitted.

Output (tab separated):
    host_chr  host_pos  host_strand  hbv_pos  hbv_strand  nss  read_id
    host_mapq  hbv_mapq  event_id  event_type

Conventions:
    - coordinates are 1-based, on the + strand of each reference
    - host_pos/hbv_pos are the genomic positions of the last host base and
      the first HBV base adjacent to the junction
    - secondary host_pos is the 1-based start of the alternative alignment
    - event_id is a unique per-file id (<sample>_bpNNN); event_type is
      'primary' or 'secondary'
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bwa_util

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def parse_cigar(c):
    return [(int(n), op) for n, op in CIGAR_RE.findall(c)]


def aln_interval(rname, pos, cigar, strand):
    """q_start/q_end: 0-based, exclusive query coordinates.
    r_start: 1-based; r_end_incl: 1-based, inclusive.
    Trailing soft/hard clips do NOT extend q_end."""
    q = 0
    ref = pos - 1
    q_start = None
    q_end = None
    for n, op in parse_cigar(cigar):
        if op in 'SH':
            q += n
        elif op in 'MI=X':
            if q_start is None:
                q_start = q
            q += n
            ref += n
            q_end = q
        elif op == 'I':
            q += n
            if q_start is not None:
                q_end = q
        elif op in 'DN':
            ref += n
    ref_len = ref - (pos - 1)
    return {
        'rname': rname,
        'strand': strand,
        'q_start': q_start,
        'q_end': q_end if q_end is not None else q,
        'r_start': pos,
        'r_end_incl': pos + ref_len - 1,
        'ref_len': ref_len,
    }


def ref_pos(a, q):
    if a['strand'] == '+':
        return a['r_start'] + (q - a['q_start'])
    return a['r_end_incl'] - (q - a['q_start'])


def ridx_of(flag):
    # bwa mem strips the /1 /2 suffix from read names, so distinguish the
    # two mates by the SAM pair flag instead:
    #   0x40 (64) = first in pair, 0x80 (128) = second in pair
    if flag & 0x40:
        return 0
    if flag & 0x80:
        return 1
    return 0


def collect(path):
    """Return {base: {'len': {ridx: read_len}, 'alns': [(ridx, aln)]}}."""
    reads = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            qname = f[0]
            base = qname[:-2] if (qname.endswith('/1') or qname.endswith('/2')) else qname
            flag = int(f[1])
            if flag & 0x800:  # supplementary alignment (SAM flag 2048) -> ignore
                continue
            rname = f[2]
            pos = int(f[3])
            cigar = f[5]
            if rname == '*' or cigar == '*':
                continue
            strand = '-' if (flag & 0x10) else '+'
            seq = f[9] if len(f) > 9 else ''
            read_len = len(seq)
            mapq = int(f[4]) if len(f) > 4 else 0
            ridx = ridx_of(flag)
            rd = reads.setdefault(base, {'len': {}, 'alns': [], 'seq': {}})
            rd['len'][ridx] = max(rd['len'].get(ridx, 0), read_len)
            if seq:
                rd['seq'][ridx] = seq
            aln = aln_interval(rname, pos, cigar, strand)
            aln['mapq'] = mapq
            rd['alns'].append((ridx, aln))
    return reads


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host_sam')
    p.add_argument('hbv_sam')
    p.add_argument('--min-match', type=int, default=30,
                   help='minimum match length (bp) on host and HBV')
    p.add_argument('--host-index', default='',
                   help='host bwa index prefix; when given, reads whose host '
                        'alignment has MAPQ == 0 are re-mapped with bwa mem -a '
                        'to report ALL possible host positions (secondary events)')
    p.add_argument('--min-score', type=int, default=15,
                   help='bwa mem -T score threshold for the secondary re-map')
    p.add_argument('--sample', default='sample',
                   help='sample id used as the event_id prefix')
    args = p.parse_args()

    host = collect(args.host_sam)
    hbv = collect(args.hbv_sam)

    primaries = []
    for base in host:
        if base not in hbv:
            continue
        hrecs = host[base]
        vrecs = hbv[base]

        # common-coordinate offsets so the two ends of a pair share one axis
        lens = dict(hrecs['len'])
        for k, v in vrecs['len'].items():
            lens[k] = max(lens.get(k, 0), v)
        off = {}
        total = 0
        for ridx in sorted(lens):
            off[ridx] = total
            total += lens[ridx]

        host_alns = hrecs['alns']
        hbv_alns = vrecs['alns']
        best_host = max(host_alns, key=lambda t: t[1]['ref_len'])
        best_hbv = max(hbv_alns, key=lambda t: t[1]['ref_len'])
        if best_host[1]['ref_len'] < args.min_match or best_hbv[1]['ref_len'] < args.min_match:
            continue

        host_lo = min(a['q_start'] + off[r] for r, a in host_alns)
        host_hi = max(a['q_end'] + off[r] for r, a in host_alns)
        hbv_lo = min(a['q_start'] + off[r] for r, a in hbv_alns)
        hbv_hi = max(a['q_end'] + off[r] for r, a in hbv_alns)

        def aln_at(alns, q):
            cands = [(r, a) for (r, a) in alns
                     if a['q_start'] + off[r] <= q < a['q_end'] + off[r]]
            if not cands:
                return None
            return max(cands, key=lambda t: t[1]['ref_len'])

        if host_hi <= hbv_lo:
            # host upstream -> HBV downstream
            q_host = host_hi - 1   # last host base
            q_hbv = hbv_lo         # first HBV base
        elif hbv_hi <= host_lo:
            # HBV upstream -> host downstream
            q_hbv = hbv_hi - 1
            q_host = host_lo
        else:
            # host and HBV alignments overlap on the read (assembled junction
            # reads usually share a few bases at the junction): use the overlap
            # midpoint as the junction position instead of dropping the read
            q_mid = (max(host_lo, hbv_lo) + min(host_hi, hbv_hi)) // 2
            q_host = q_mid
            q_hbv = q_mid

        ha = aln_at(host_alns, q_host)
        va = aln_at(hbv_alns, q_hbv)
        if ha is None or va is None:
            continue
        host_pos = ref_pos(ha[1], q_host - off[ha[0]])
        hbv_pos = ref_pos(va[1], q_hbv - off[va[0]])
        primaries.append({
            'read_id': base,
            'host_seq': hrecs['seq'].get(ha[0], ''),
            'host_chr': ha[1]['rname'],
            'host_pos': host_pos,
            'host_strand': ha[1]['strand'],
            'hbv_pos': hbv_pos,
            'hbv_strand': va[1]['strand'],
            'host_mapq': ha[1]['mapq'],
            'hbv_mapq': va[1]['mapq'],
        })

    # secondary events: enumerate ALL host positions of mapq==0 reads
    alt = {}
    if args.host_index:
        mapq0 = [(pr['read_id'], pr['host_seq']) for pr in primaries
                 if pr['host_mapq'] == 0 and pr['host_seq']]
        alt = bwa_util.remap_host_all(mapq0, args.host_index, args.min_score)

    out = sys.stdout
    out.write('host_chr\thost_pos\thost_strand\thbv_pos\thbv_strand\tnss\tread_id\thost_mapq\thbv_mapq\tevent_id\tevent_type\n')
    n = 0
    for pr in primaries:
        n += 1
        eid = '%s_bp%03d' % (args.sample, n)  # one id per PRIMARY event
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t1\t{5}\t{6}\t{7}\t{8}\tprimary\n'.format(
            pr['host_chr'], pr['host_pos'], pr['host_strand'],
            pr['hbv_pos'], pr['hbv_strand'], pr['read_id'],
            pr['host_mapq'], pr['hbv_mapq'], eid))
        if pr['host_mapq'] == 0:
            # secondary events share the SAME event_id as their primary
            for (is_sec, rname, pos, strand, mapq, rlen) in alt.get(pr['read_id'], []):
                if not is_sec or rlen < args.min_match:
                    continue
                if rname == pr['host_chr'] and abs(pos - pr['host_pos']) <= 5:
                    continue  # same locus as the primary -> skip
                out.write('{0}\t{1}\t{2}\t{3}\t{4}\t1\t{5}\t{6}\t{7}\t{8}\tsecondary\n'.format(
                    rname, pos, strand, pr['hbv_pos'], pr['hbv_strand'],
                    pr['read_id'], mapq, pr['hbv_mapq'], eid))


if __name__ == '__main__':
    main()
