#!/usr/bin/env python3
"""Prefer sequenced junctions for discordant support; otherwise retain estimates.

Both human and HBV distances and human-to-HBV traversal orientations must agree.
--dist bounds the total unsequenced span (human distance + viral distance).
"""
import argparse
import sys
from junctions import *
from call_breakpoints import call_rows
from merge_breakpoints import cluster


def attach_pairs(sequenced, pairs, dist):
    unresolved=[]
    for pair in pairs:
        placements=[pair]+pair['secondary_events']
        matches=[]
        for target in sequenced:
            supporting=[a for a in placements if pair_supports(a,target,dist)]
            if supporting:
                best=min(supporting,key=lambda a:abs(a['host_pos']-target['host_pos'])+abs(a['hbv_pos']-target['hbv_pos']))
                matches.append((abs(best['host_pos']-target['host_pos'])+abs(best['hbv_pos']-target['hbv_pos']),target))
        if not matches:
            unresolved.append(pair)
            continue
        matches.sort(key=lambda x:(x[0],key(x[1])))
        target=matches[0][1]
        if len(matches)>1:
            for evidence in pair['molecules'].values():
                evidence['ambiguous']=True
        combine_molecules(target['molecules'],pair['molecules'])
        # Keep alternatives attached to the supporting molecule, without
        # replacing a sequenced junction with the estimated mate-end position.
        for mid in pair['molecules']:
            target['molecules'][mid]['alternative_placements']=pair['secondary_events']
            if len(matches)>1:
                target['molecules'][mid]['compatible_events']=[t['event_id'] for _,t in matches]
        refresh(target)
    return sequenced+cluster(unresolved,dist)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('host_sam');p.add_argument('hbv_sam');p.add_argument('breakpoints')
    p.add_argument('--dist',type=int,default=300)
    p.add_argument('--merge-dist',type=int,default=20)
    add_arguments(p)
    args=p.parse_args()
    with open(args.breakpoints) as handle:
        merged=read_rows(handle)
    unmerged=call_rows(collect(args.host_sam),collect(args.hbv_sam),args,'unmerged',True)
    sequenced=cluster(merged+[r for r in unmerged if r['breakpoint_resolution']!='estimated_pair'],args.merge_dist)
    pairs=[r for r in unmerged if r['breakpoint_resolution']=='estimated_pair']
    write_rows(attach_pairs(sequenced,pairs,args.dist),sys.stdout)

if __name__=='__main__':
    main()
