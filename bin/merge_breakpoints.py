#!/usr/bin/env python3
"""Cluster compatible junctions on BOTH references; retain molecule ambiguity."""
import argparse
import sys
from junctions import *


def merge_into(target, row):
    combine_molecules(target['molecules'],row['molecules'])
    alternatives = {key(a): a for a in target['secondary_events']+row['secondary_events']}
    alternatives.pop(key(target),None)
    target['secondary_events'] = [alternatives[k] for k in sorted(alternatives)]
    for prefix in ('host','hbv'):
        target[prefix+'_ci_start'] = min(target[prefix+'_ci_start'],row[prefix+'_ci_start'])
        target[prefix+'_ci_end'] = max(target[prefix+'_ci_end'],row[prefix+'_ci_end'])
    if row['source']=='merged':
        target['source']='merged'
    refresh(target)


def cluster(rows, dist):
    # Fixed representative and complete-link bounds prevent transitive chains
    # joining sites farther apart than the requested distance.
    groups=[]
    ordered=sorted(rows,key=lambda r:(-r['unique_nss'],-r['nss'],-r['host_mapq']-r['hbv_mapq'],key(r),r['event_id']))
    for row in ordered:
        group=next((g for g in groups if all(compatible(row,x,dist) and
                    (row['breakpoint_resolution']=='estimated_pair') == (x['breakpoint_resolution']=='estimated_pair')
                    for x in g)),None)
        if group is None:
            groups.append([row])
        else:
            group.append(row)
    result=[]
    for group in groups:
        rep=group[0]
        for row in group[1:]:
            merge_into(rep,row)
        result.append(rep)
    return sorted(result,key=key)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dist',type=int,default=20)
    args=p.parse_args()
    write_rows(cluster(read_rows(sys.stdin),args.dist),sys.stdout)

if __name__=='__main__':
    main()
