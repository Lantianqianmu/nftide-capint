#!/usr/bin/env python3
"""Call sequenced junctions using shared strand/CIGAR-aware geometry."""
import argparse
import sys
from types import SimpleNamespace
import bwa_util
from junctions import *


def augment_host(host, args):
    # Remap full original sequences, not reference-oriented SAM sequences.
    queries, lookup = [], {}
    for name, alns in host.items():
        for mate in sorted({a.mate for a in alns}):
            subset = [a for a in alns if a.mate == mate]
            if not any(a.mapq < args.min_mapq for a in subset):
                continue
            full = next((a for a in subset if len(a.sequence) == a.qlen), None)
            if full:
                token = 'alt%d' % len(lookup)
                lookup[token] = (name,mate)
                queries.append((token,full.sequence))
    for token, alns in bwa_util.remap_host_all(queries,args.host_index,args.min_score).items():
        name,mate = lookup[token]
        for a in alns:
            a.name, a.mate = name,mate
            host[name].append(a)


def candidates(hs, vs, args, allow_pairs=False):
    found = []
    for h in hs:
        for v in vs:
            options = SimpleNamespace(**vars(args))
            options.best_single_score = max((a.score for a in hs+vs if a.mate == h.mate),default=0)
            row = split_junction(h,v,options)
            if row:
                found.append((h.score+v.score,row))
    if not found and allow_pairs:
        # A discordant pair needs different mates, almost full-length anchors,
        # and no comparable alignment of either mate to the other genome.
        for h in hs:
            for v in vs:
                if h.mate == v.mate or min(h.matched,v.matched) <= args.min_match:
                    continue
                if h.matched/h.qlen < .9 or v.matched/v.qlen < .9:
                    continue
                if any(a.mate == h.mate and a.score >= h.score-args.score_margin for a in vs):
                    continue
                if any(a.mate == v.mate and a.score >= v.score-args.score_margin for a in hs):
                    continue
                hp,vp = h.boundary(),v.boundary()
                found.append((h.score+v.score,location(h,v,hp,vp,h.strand,flip(v.strand),'estimated_pair')))
    if not found:
        return []
    best = max(score for score,_ in found)
    unique = {}
    for score,row in sorted(found,key=lambda x: (-x[0],-x[1]['host_mapq'],-x[1]['hbv_mapq'],key(x[1]))):
        if score >= best-args.alt_score_delta:
            unique.setdefault(key(row),row)
    return list(unique.values())


def call_rows(host, virus, args, source, allow_pairs=False):
    augment_host(host,args)
    rows = []
    for name in sorted(host.keys() & virus.keys()):
        hs,vs = host[name],virus[name]
        options = candidates(hs,vs,args,allow_pairs)
        if not options:
            continue
        row = options[0]
        ambiguous = len(options)>1 or min(row['host_mapq'],row['hbv_mapq']) < args.min_mapq
        seqs = []
        for mate in sorted({a.mate for a in hs+vs}):
            full = next((a.sequence for a in hs+vs if a.mate==mate and len(a.sequence)==a.qlen), '')
            seqs.append(full)
        mid = molecule_id(name,seqs)
        row.update(source=source,event_id=args.sample+'_'+mid,event_type='primary',
                   secondary_events=options[1:],molecules={mid:dict(ambiguous=ambiguous,
                   kind='pair' if row['breakpoint_resolution']=='estimated_pair' else 'split',
                   copies=raw_count(name),junction_gap_length=row['junction_gap_length'])})
        rows.append(refresh(row))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host_sam'); p.add_argument('hbv_sam')
    add_arguments(p)
    args=p.parse_args()
    write_rows(call_rows(collect(args.host_sam),collect(args.hbv_sam),args,'merged'),sys.stdout)

if __name__ == '__main__':
    main()
