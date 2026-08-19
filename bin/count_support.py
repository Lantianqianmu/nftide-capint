#!/usr/bin/env python3
"""Combine merged-read breakpoints with unmerged chimeric-pair support.

Unmerged (not spliced by bbmerge) chimeric pairs do not contain a breakpoint
by themselves, but a pair whose one end maps to the host and the other to HBV
flanks an integration junction. Such pairs:

  (a) add NSS support to a nearby breakpoint already defined by merged reads, and
  (b) define a breakpoint of their own when no merged-read breakpoint is nearby
      (clustered within --dist, so several pairs at the same locus merge into
      one breakpoint with NSS = number of pairs).

A pair is a valid support as soon as it has a host match and an HBV match of
>= --min-match bp each (host and HBV taken independently as the longest
primary alignment of the pair to each genome).

Two situations are handled, and BOTH report an EXACT breakpoint coordinate
(not just the start of the longest alignment):

  1. CHIMERIC MATE - one mate contains both a host and an HBV segment
     (e.g. host CIGAR "64S85M" + HBV CIGAR "81S68M" on the same read). The
     junction is resolved inside that read from the two CIGARs and the exact
     host / HBV coordinates at the junction are reported.

  2. SPLIT PAIR - one mate is pure host and the other pure HBV. The junction
     lies between the two reads, i.e. at the 3' (fragment-interior) border of
     each mate; those border coordinates are reported.

If the exact coordinate falls within --dist of a breakpoint already defined by
merged reads, the pair adds NSS support to it instead of creating a new row.

The output is a COMBINED breakpoint table (same columns as MERGE_SIGNAL plus a
source column and the MAPQ of the representative read/pair):
    host_chr  host_pos  host_strand  hbv_pos  hbv_strand  nss
    source  host_mapq  hbv_mapq

Inputs:
    host.sam / hbv.sam   BWA-MEM alignments of the unmerged chimeric pairs
                         against the host and HBV genomes
    breakpoints          raw breakpoint table (from REMAP; one row per
                         merged/spliced read)

Only primary alignments are used; supplementary alignments (SAM flag 2048)
are ignored.
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


class Aln(object):
    """Primary alignment of one mate, normalized to forward-read coordinates.

    q_start / q_end : 0-based half-open query interval in the FORWARD read
                      (reverse-strand alignments are flipped so the host and
                      HBV segments of the same mate share one coordinate axis).
    r_start         : 1-based reference start (leftmost matched base)
    r_end_incl      : 1-based inclusive reference end
    ref_len         : reference bases covered (M/=/X + D)
    """
    __slots__ = ('rname', 'r_start', 'r_end_incl', 'strand', 'q_start', 'q_end',
                 'ref_len', 'mapq', 'ridx')

    def __init__(self, rname, r_start, r_end_incl, strand, q_start, q_end,
                 ref_len, mapq, ridx):
        self.rname = rname
        self.r_start = r_start
        self.r_end_incl = r_end_incl
        self.strand = strand
        self.q_start = q_start
        self.q_end = q_end
        self.ref_len = ref_len
        self.mapq = mapq
        self.ridx = ridx

    def ref_pos(self, q):
        """1-based reference position of forward-read position q (0-based)."""
        if self.strand == '+':
            return self.r_start + (q - self.q_start)
        return self.r_end_incl - (q - self.q_start)

    def border(self):
        """3' border of this mate - the end facing the fragment interior."""
        return self.r_end_incl if self.strand == '+' else self.r_start


def exact_junction(h, v):
    """Exact breakpoint of a chimeric mate.

    h / v : Aln of the host / HBV segment on the SAME mate.
    Returns (host_pos, hbv_pos): 1-based reference coordinates of the two
    junction bases (last host base / first HBV base adjacent to the junction;
    when the two segments overlap on the read, the overlap midpoint is used).
    """
    h_lo, h_hi = h.q_start, h.q_end
    v_lo, v_hi = v.q_start, v.q_end
    if h_hi <= v_lo:
        # host upstream -> HBV downstream: junction at host 3' / HBV 5'
        q_host = h_hi - 1
        q_hbv = v_lo
    elif v_hi <= h_lo:
        # HBV upstream -> host downstream: junction at HBV 3' / host 5'
        q_hbv = v_hi - 1
        q_host = h_lo
    else:
        # segments overlap on the read (junction reads usually share bases)
        q_mid = (max(h_lo, v_lo) + min(h_hi, v_hi)) // 2
        q_host = q_mid
        q_hbv = q_mid
    return h.ref_pos(q_host), v.ref_pos(q_hbv)


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
    """Return {base: {'len': {ridx: read_len}, 'alns': [(ridx, Aln)]}}.

    Query intervals are normalized to forward-read coordinates so the host and
    HBV segments of the same mate are directly comparable.
    """
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

            # query interval in the alignment's own (possibly rev-comp) space
            q = 0
            ref = pos - 1
            qs = qe = None
            for n, op in parse_cigar(cigar):
                if op in 'SH':
                    q += n
                elif op in 'MI=X':
                    if qs is None:
                        qs = q
                    q += n
                    ref += n
                    qe = q
                elif op == 'I':
                    q += n
                elif op in 'DN':
                    ref += n
            if qs is None:  # no aligned base (all soft/hard clipped)
                continue
            ref_len = ref - (pos - 1)
            # normalize to forward-read coordinates
            if strand == '-':
                q_start = read_len - qe
                q_end = read_len - qs
            else:
                q_start = qs
                q_end = qe

            aln = Aln(rname, pos, pos + ref_len - 1, strand, q_start, q_end,
                      ref_len, mapq, ridx)
            rd = reads.setdefault(base, {'len': {}, 'alns': [], 'seq': {}})
            rd['len'][ridx] = max(rd['len'].get(ridx, 0), read_len)
            if seq:
                rd['seq'][ridx] = seq
            rd['alns'].append((ridx, aln))
    return reads


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host_sam')
    p.add_argument('hbv_sam')
    p.add_argument('breakpoints')
    p.add_argument('--min-match', type=int, default=30,
                   help='minimum match length (bp) on host and HBV')
    p.add_argument('--dist', type=int, default=300,
                   help='max distance (bp) for clustering unmerged pairs and '
                        'matching them to merged breakpoints (host side)')
    p.add_argument('--host-index', default='',
                   help='host bwa index prefix; when given, unmerged-defined '
                        'breakpoints with host MAPQ == 0 are re-mapped with '
                        'bwa mem -a to report ALL possible host positions '
                        '(secondary events)')
    p.add_argument('--min-score', type=int, default=15,
                   help='bwa mem -T score threshold for the secondary re-map')
    p.add_argument('--sample', default='sample',
                   help='sample id used as the event_id prefix')
    args = p.parse_args()

    # load breakpoints (raw REMAP output: ... nss read_id host_mapq hbv_mapq)
    merged = []  # each: [chr, pos, hstrand, hbv_pos, vstrand, nss, nbp, host_mapq, hbv_mapq]
    with open(args.breakpoints) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            toks = line.split('\t')
            if len(toks) < 7:
                continue
            try:
                mp = int(toks[1])
                mv = int(toks[3])
                ms = int(toks[5])
            except ValueError:
                continue  # skip header / non-numeric rows
            # col 6 is n_breakpoints (merged table) or read_id (raw table)
            try:
                mn = int(toks[6])
            except ValueError:
                mn = 1  # raw per-read row: one read == one support
            merged.append([toks[0], mp, toks[2], mv, toks[4], ms, mn,
                           int(toks[7]) if len(toks) > 7 and toks[7] else 0,
                           int(toks[8]) if len(toks) > 8 and toks[8] else 0,
                           toks[9] if len(toks) > 9 else '',
                           toks[10] if len(toks) > 10 else 'primary'])

    host = collect(args.host_sam)
    hbv = collect(args.hbv_sam)

    # one breakpoint per chimeric unmerged pair (EXACT coordinate):
    #   chimeric mate -> exact junction resolved from the two CIGARs on that mate
    #   split pair    -> 3' (fragment-interior) border of the host and HBV mate
    pairs = []  # (chr, host_pos, hstrand, hbv_pos, vstrand, host_mapq, hbv_mapq)
    for base in host:
        if base not in hbv:
            continue
        hrecs = host[base]
        vrecs = hbv[base]

        # longest primary host & HBV alignment of each mate
        by_ridx = {}
        for ridx, a in hrecs['alns']:
            m = by_ridx.setdefault(ridx, {'host': None, 'hbv': None})
            if m['host'] is None or a.ref_len > m['host'].ref_len:
                m['host'] = a
        for ridx, a in vrecs['alns']:
            m = by_ridx.setdefault(ridx, {'host': None, 'hbv': None})
            if m['hbv'] is None or a.ref_len > m['hbv'].ref_len:
                m['hbv'] = a

        mates = []
        for ridx, m in by_ridx.items():
            h = m['host'] if (m['host'] and m['host'].ref_len >= args.min_match) else None
            v = m['hbv'] if (m['hbv'] and m['hbv'].ref_len >= args.min_match) else None
            if h is not None or v is not None:
                mates.append((h, v))

        # (1) chimeric mate: both segments on the same read -> exact junction
        chim = [m for m in mates if m[0] is not None and m[1] is not None]
        if chim:
            chim.sort(key=lambda m: -m[0].ref_len)  # mate with the longest host
            h, v = chim[0]
            hpos, vpos = exact_junction(h, v)
            pairs.append((h.rname, hpos, h.strand, vpos, v.strand,
                          h.mapq, v.mapq, base,
                          hrecs['seq'].get(h.ridx, '')))
            continue

        # (2) split pair: one mate pure host, the other pure HBV -> 3' border
        host_only = [m[0] for m in mates if m[0] is not None]
        hbv_only = [m[1] for m in mates if m[1] is not None]
        if host_only and hbv_only:
            h = max(host_only, key=lambda a: a.ref_len)
            v = max(hbv_only, key=lambda a: a.ref_len)
            hpos = h.border()
            vpos = v.border()
            pairs.append((h.rname, hpos, h.strand, vpos, v.strand,
                          h.mapq, v.mapq, base,
                          hrecs['seq'].get(h.ridx, '')))
            continue
        # otherwise: no valid host+HBV combination on this read

    # cluster unmerged pairs by host locus (within --dist)
    pairs.sort(key=lambda x: (x[0], x[1]))
    clusters = []
    for pr in pairs:
        if clusters:
            cl = clusters[-1]
            if pr[0] == cl[0] and abs(pr[1] - cl[1]) <= args.dist:
                cl[2].append(pr)
                continue
        clusters.append([pr[0], pr[1], [pr]])

    # add unmerged support to merged breakpoints; else make a new breakpoint
    nss = [int(m[5]) for m in merged]
    new_bps = []  # unmerged-defined breakpoints

    for (hchr, rep_pos, plist) in clusters:
        # nearest merged breakpoint on the same chromosome within --dist
        best_i = -1
        best_d = None
        for i, m in enumerate(merged):
            if m[0] != hchr:
                continue
            d = abs(m[1] - rep_pos)
            if d <= args.dist and (best_d is None or d < best_d):
                best_d = d
                best_i = i
        if best_i >= 0:
            nss[best_i] += len(plist)
        else:
            rep = plist[len(plist) // 2]  # representative pair of the cluster
            new_bps.append([hchr, rep_pos, rep[2], rep[3], rep[4],
                            len(plist), len(plist), rep[5], rep[6],
                            rep[7], rep[8]])

    # event ids: continue after the highest id seen in the raw (merged) rows
    eid_n = 0
    for m in merged:
        mm = re.search(r'(\d+)$', m[9] or '')
        if mm:
            eid_n = max(eid_n, int(mm.group(1)))

    # secondary events for unmerged-defined breakpoints with host_mapq == 0
    alt = {}
    if args.host_index:
        mapq0 = [(nb[9], nb[10]) for nb in new_bps
                 if nb[7] == 0 and nb[10]]
        alt = bwa_util.remap_host_all(mapq0, args.host_index, args.min_score)

    out = sys.stdout
    out.write('host_chr\thost_pos\thost_strand\thbv_pos\thbv_strand\tnss\tsource\thost_mapq\thbv_mapq\tevent_id\tevent_type\n')
    for i, m in enumerate(merged):
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tmerged\t{6}\t{7}\t{8}\t{9}\n'.format(
            m[0], m[1], m[2], m[3], m[4], nss[i], m[7], m[8],
            m[9], m[10]))
    for nb in new_bps:
        eid_n += 1
        eid = '%s_bp%03d' % (args.sample, eid_n)  # one id per PRIMARY event
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tunmerged\t{6}\t{7}\t{8}\tprimary\n'.format(
            nb[0], nb[1], nb[2], nb[3], nb[4], nb[5], nb[7], nb[8], eid))
        if nb[7] == 0:
            # secondary events share the SAME event_id as their primary
            # req(2): a secondary of an unassembled event that falls within
            # --dist of a PE-assembled (merged) breakpoint merges into it
            # instead of creating a separate row. At event-id level the SAME
            # read/event counts at most once per merged breakpoint it overlaps.
            absorbed = set()
            for (is_sec, rname, pos, strand, mapq, rlen) in alt.get(nb[9], []):
                if not is_sec or rlen < args.min_match:
                    continue
                if rname == nb[0] and abs(pos - nb[1]) <= 5:
                    continue
                j = -1
                jd = None
                for i, m in enumerate(merged):
                    if m[0] != rname:
                        continue
                    d = abs(m[1] - pos)
                    if d <= args.dist and (jd is None or d < jd):
                        jd = d
                        j = i
                if j >= 0:
                    absorbed.add(j)   # +1 NSS per distinct merged breakpoint
                    continue
                out.write('{0}\t{1}\t{2}\t{3}\t{4}\t1\tunmerged\t{5}\t{6}\t{7}\tsecondary\n'.format(
                    rname, pos, strand, nb[3], nb[4], mapq, nb[8], eid))
            for j in absorbed:
                nss[j] += 1


if __name__ == '__main__':
    main()
