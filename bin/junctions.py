#!/usr/bin/env python3
"""Shared SAM geometry and evidence tables. Coordinates are 1-based.

Strands describe traversal from human into HBV, independent of read direction.
Molecule IDs are sequence-deduplication proxies, not UMI measurements.
"""
import csv
from table_io import dict_reader
import hashlib
import json
import re
from dataclasses import dataclass

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
FIELDS = 'host_chr host_pos host_strand hbv_pos hbv_strand nss source host_mapq hbv_mapq event_id event_type breakpoint_resolution host_ci_start host_ci_end hbv_ci_start hbv_ci_end unique_nss ambiguous_nss split_nss pair_nss raw_read_pairs secondary_events molecules junction_gap_length'.split()

def rc(s):
    return s.translate(COMP)[::-1]

def flip(s):
    return '-' if s == '+' else '+'

def molecule_id(name, seqs):
    match = re.search(r'\|mol=([a-f0-9]+)', name)
    if match:
        return match[1]
    seqs = tuple(s.upper() for s in seqs)
    reverse = (rc(seqs[0]),) if len(seqs) == 1 else tuple(reversed(seqs))
    return hashlib.sha256('|'.join(min(seqs, reverse)).encode()).hexdigest()[:32]

def raw_count(name):
    match = re.search(r'\|copies=(\d+)', name)
    return int(match[1]) if match else 1

@dataclass
class Alignment:
    name: str
    mate: int
    chrom: str
    strand: str
    mapq: int
    score: int
    secondary: bool
    qlen: int
    positions: dict
    sequence: str

    @property
    def lo(self):
        return min(self.positions)

    @property
    def hi(self):
        return max(self.positions) + 1

    @property
    def matched(self):
        return len(self.positions)

    def boundary(self, inward=True):
        return self.positions[max(self.positions) if inward else min(self.positions)]

def parse_alignment(f):
    flag = int(f[1])
    if flag & 4 or f[2] == '*' or f[5] == '*':
        return None
    ops = [(int(n), op) for n, op in CIGAR.findall(f[5])]
    length = sum(n for n, op in ops if op in 'MISH=X')
    q, ref, points = 0, int(f[3]), {}
    for n, op in ops:
        if op in 'M=X':
            for i in range(n):
                points[length - 1 - q - i if flag & 16 else q + i] = ref + i
            q += n
            ref += n
        elif op in 'ISH':
            q += n
        elif op in 'DN':
            ref += n
    if not points:
        return None
    tags = {t[:2]: t[5:] for t in f[11:]}
    seq = '' if f[9] == '*' else f[9]
    if flag & 16:
        seq = rc(seq)
    return Alignment(f[0], 1 if flag & 128 else 0, f[2], '-' if flag & 16 else '+',
                     int(f[4]), int(tags.get('AS', len(points))), bool(flag & 256),
                     length, points, seq)

def collect(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            if len(f) < 11:
                continue
            aln = parse_alignment(f)
            if aln:
                result.setdefault(aln.name, []).append(aln)
    return result

def location(h, v, hp, vp, hs, vs, resolution, hci=None, vci=None, junction_gap_length="NA"):
    return dict(host_chr=h.chrom, host_pos=hp, host_strand=hs, hbv_pos=vp,
                hbv_strand=vs, host_mapq=h.mapq, hbv_mapq=v.mapq,
                breakpoint_resolution=resolution, junction_gap_length=junction_gap_length,
                host_ci_start=(hci or (hp, hp))[0], host_ci_end=(hci or (hp, hp))[1],
                hbv_ci_start=(vci or (vp, vp))[0], hbv_ci_end=(vci or (vp, vp))[1])

def split_junction(h, v, args):
    if h.mate != v.mate or h.qlen != v.qlen:
        return None
    if min(h.matched, v.matched) <= args.min_match:
        return None
    # Each genome must explain sequence not covered by the other alignment.
    if min(len(set(h.positions) - set(v.positions)),
           len(set(v.positions) - set(h.positions))) < args.min_unique:
        return None
    left, right = (h, v) if h.lo < v.lo else (v, h)
    if left.hi >= right.hi:
        return None
    gap = right.lo - left.hi
    if gap < -args.max_overlap:
        return None
    # Overlapping explanations require additional score evidence. Disjoint
    # anchors longer than min_match are accepted regardless of gap length.
    if gap < 0 and h.score + v.score + gap < args.best_single_score + args.score_margin:
        return None
    if gap < 0:
        cut = (right.lo + left.hi) // 2
        lq = max((q for q in left.positions if q < cut), default=None)
        rq = min((q for q in right.positions if q >= cut), default=None)
    else:
        lq, rq = left.hi - 1, right.lo
    if lq is None or rq is None:
        return None
    lp, rp = left.positions[lq], right.positions[rq]
    hp, vp = (lp, rp) if left is h else (rp, lp)
    hs, vs = (h.strand, v.strand) if left is h else (flip(h.strand), flip(v.strand))
    resolution = 'sequenced' if gap == 0 else ('microhomology' if gap < 0 else 'sequenced_gap')
    hci, vci = (hp, hp), (vp, vp)
    if gap < 0:
        # Include both flanking choices in the microhomology uncertainty range.
        def span(a, pos):
            vals = [r for q, r in a.positions.items() if right.lo - 1 <= q <= left.hi]
            return min(vals + [pos]), max(vals + [pos])
        hci, vci = span(h, hp), span(v, vp)
    return location(h, v, hp, vp, hs, vs, resolution, hci, vci, junction_gap_length=max(0, gap))

def key(row):
    return tuple(row[k] for k in ('host_chr', 'host_pos', 'host_strand', 'hbv_pos', 'hbv_strand'))

def orientation(row):
    return row['host_chr'], row['host_strand'], row['hbv_strand']

def compatible(a, b, dist):
    return orientation(a) == orientation(b) and all(abs(int(a[k])-int(b[k])) <= dist for k in ('host_pos','hbv_pos'))

def pair_supports(pair, split, dist, slop=10):
    if orientation(pair) != orientation(split):
        return False
    hd = (int(split['host_pos'])-int(pair['host_pos'])) * (1 if pair['host_strand']=='+' else -1)
    vd = (int(pair['hbv_pos'])-int(split['hbv_pos'])) * (1 if pair['hbv_strand']=='+' else -1)
    return -slop <= hd <= dist and -slop <= vd <= dist and hd + vd <= dist

def add_arguments(p):
    p.add_argument('--min-match', type=int, default=30,
                   help='each anchor must have strictly more than this many aligned bases')
    p.add_argument('--min-unique', type=int, default=20)
    p.add_argument('--max-overlap', type=int, default=10)
    p.add_argument('--score-margin', type=int, default=10)
    p.add_argument('--min-mapq', type=int, default=20)
    p.add_argument('--alt-score-delta', type=int, default=5)
    p.add_argument('--host-index', default='')
    p.add_argument('--min-score', type=int, default=15)
    p.add_argument('--sample', default='sample')

def refresh(row):
    mols = row['molecules']
    row['nss'] = len(mols)
    row['unique_nss'] = sum(not m['ambiguous'] for m in mols.values())
    row['ambiguous_nss'] = row['nss'] - row['unique_nss']
    row['split_nss'] = sum(m['kind'] == 'split' for m in mols.values())
    row['pair_nss'] = row['nss'] - row['split_nss']
    row['raw_read_pairs'] = sum(m['copies'] for m in mols.values())
    row['event_type'] = 'ambiguous' if not row['unique_nss'] else 'primary'
    return row

def combine_molecules(target, source):
    for mid, evidence in source.items():
        if mid not in target:
            target[mid] = dict(evidence)
        else:
            old = target[mid]
            old['copies'] = max(old['copies'], evidence['copies'])
            old['ambiguous'] = old['ambiguous'] or evidence['ambiguous']
            if evidence['kind'] == 'split':
                if old['kind'] == 'pair':
                    old['junction_gap_length'] = evidence.get('junction_gap_length', 'NA')
                old['kind'] = 'split'

def read_rows(handle):
    rows = []
    for row in dict_reader(handle):
        if 'molecules' not in row or 'junction_gap_length' not in row:
            raise ValueError('Old breakpoint schema: rerun REMAP and downstream stages')
        for k in ('host_pos','hbv_pos','host_mapq','hbv_mapq','host_ci_start','host_ci_end','hbv_ci_start','hbv_ci_end'):
            row[k] = int(row[k])
        for k in ('molecules','secondary_events'):
            row[k] = json.loads(row[k])
        if row['junction_gap_length'] != 'NA':
            row['junction_gap_length'] = int(row['junction_gap_length'])
        rows.append(refresh(row))
    return rows

def write_rows(rows, handle):
    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    for row in rows:
        out = dict(refresh(row))
        for k in ('molecules','secondary_events'):
            out[k] = json.dumps(out[k], sort_keys=True, separators=(',', ':'))
        writer.writerow(out)
