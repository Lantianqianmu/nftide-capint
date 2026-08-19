#!/usr/bin/env python3
"""Merge HBV integration breakpoints within DIST bp (HIVID step 5).

Runs on the COMBINED breakpoint table (after unmerged support counting and
secondary-event calling). Breakpoints whose host position lies within DIST bp
are merged into one event:

  - if any member of a cluster is a PRIMARY event, the merged event is primary
    and the member with the highest NSS among the primaries is kept (its
    coordinates, MAPQ and event_id); secondary events merge INTO the primary;
  - otherwise (secondary-only cluster) the member with the highest NSS is kept.

NSS is counted at the EVENT-ID level: each distinct event_id in the cluster
contributes its max NSS (so the 6bp-frame-shifted secondary hits of ONE read
count at most 1, and merging with a DIFFERENT event adds that event's support).

Event-level merge by overlapping secondary (2026-08-18): in addition to the
host-position rule above, two DIFFERENT events are merged into one event when
their HBV side (hbv_pos + hbv_strand) is identical AND their secondary loci
overlap (secondary rows land in the same host-position cluster). Overlapping
secondary hits mean the two reads may come from the same integration region
(all primary/secondary rows here have host_mapq==0, because secondary calling
only runs on mapq==0 events). The merged event keeps the highest-NSS primary's
id/coordinates and NSS is summed at event-id level.

event_id is NOT re-assigned here: it was assigned once upstream (REMAP for
merged reads, COUNT_UNMERGED_NSS for unmerged-defined events) and is piped
through unchanged. A primary event and its secondary events share one
event_id; merging keeps the id of the highest-NSS (primary) member.

Input columns (combined table from COUNT_UNMERGED_NSS):
    host_chr  host_pos  host_strand  hbv_pos  hbv_strand  nss
    source  host_mapq  hbv_mapq  event_id  event_type
Output columns:
    host_chr  host_pos  host_strand  hbv_pos  hbv_strand  nss
    source  host_mapq  hbv_mapq  event_id  event_type
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dist', type=int, default=20,
                   help='merge breakpoints within this distance (bp)')
    args = p.parse_args()

    recs = []
    for line in sys.stdin:
        line = line.rstrip('\n')
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 11:
            continue
        # host_chr host_pos hstrand hbv_pos vstrand nss source hmapq vmapq eid etype
        host_chr, host_pos, hstrand, hbv_pos, vstrand, nss, source, hmapq, vmapq, eid, etype = f[:11]
        try:
            host_pos = int(host_pos)
            hbv_pos = int(hbv_pos)
            nss = int(nss)
            hmapq = int(hmapq)
            vmapq = int(vmapq)
        except ValueError:
            continue  # skip header / non-numeric rows
        recs.append([host_chr, host_pos, hstrand, hbv_pos, vstrand, nss,
                     source, hmapq, vmapq, eid, etype])
    recs.sort(key=lambda r: (r[0], r[1]))

    # (1) existing rule: merge rows whose host positions are within DIST bp
    clusters = []
    for r in recs:
        if clusters:
            cl = clusters[-1]
            rep = max(cl, key=lambda x: x[7] + x[8])  # representative = highest MAPQ
            if r[0] == rep[0] and abs(r[1] - rep[1]) <= args.dist:
                cl.append(r)
                continue
        clusters.append([r])

    # (2) event-level merge (2026-08-18): two DIFFERENT events whose HBV side
    #     (hbv_pos + hbv_strand) is identical AND whose SECONDARY loci overlap
    #     (their secondary rows landed in the same host-position cluster) are
    #     likely from the same integration region -> the two events (primaries
    #     + all their rows) are merged into one event. Note: primary and
    #     secondary rows all have host_mapq==0 here, because secondary calling
    #     only runs for mapq==0 events.
    events = {}                       # event_id -> [(cluster_idx, row), ...]
    for ci, cl in enumerate(clusters):
        for r in cl:
            events.setdefault(r[9], []).append((ci, r))

    def hbv_side(eid):
        rows = events[eid]
        for _ci, r in rows:
            if r[10] == 'primary':
                return (r[3], r[4])
        return (rows[0][1][3], rows[0][1][4])

    def prim_cluster(eid):
        for ci, r in events[eid]:
            if r[10] == 'primary':
                return ci
        return None

    cl_sec_events = [
        set(r[9] for r in cl if r[10] == 'secondary') for cl in clusters
    ]

    m = len(clusters)
    parent_c = list(range(m))

    def findc(x):
        while parent_c[x] != x:
            parent_c[x] = parent_c[parent_c[x]]
            x = parent_c[x]
        return x

    def unionc(a, b):
        ra, rb = findc(a), findc(b)
        if ra != rb:
            parent_c[rb] = ra

    for ci in range(m):
        sec_events = sorted(cl_sec_events[ci])
        for x in range(len(sec_events)):
            for y in range(x + 1, len(sec_events)):
                ea, eb = sec_events[x], sec_events[y]
                if hbv_side(ea) != hbv_side(eb):
                    continue
                pa, pb = prim_cluster(ea), prim_cluster(eb)
                if pa is not None:
                    unionc(ci, pa)
                if pb is not None:
                    unionc(ci, pb)
                if pa is not None and pb is not None:
                    unionc(pa, pb)

    # flatten clusters after event-level merging (keep deterministic order)
    merged = {}
    for ci, cl in enumerate(clusters):
        merged.setdefault(findc(ci), []).extend(cl)
    ordered = sorted(
        merged.values(),
        key=lambda cl: min((r[0], r[1]) for r in cl))

    out = sys.stdout
    out.write('host_chr\thost_pos\thost_strand\thbv_pos\thbv_strand\tnss\tsource\thost_mapq\thbv_mapq\tevent_id\tevent_type\n')
    for cl in ordered:
        prim = [x for x in cl if x[10] == 'primary']
        if prim:
            # a primary event wins: keep the id of the primary with the
            # highest NSS (its coordinates are kept as the merged event)
            rep = max(prim, key=lambda x: (x[5], x[7] + x[8]))
        else:
            # secondary-only cluster: keep the id with the highest NSS
            rep = max(cl, key=lambda x: (x[5], x[7] + x[8]))
        # NSS at EVENT-ID level: each distinct event_id contributes its max NSS
        # (same-event rows / 6bp-frame-shift secondary hits of one read count once).
        contrib = {}
        for x in cl:
            if x[5] > contrib.get(x[9], 0):
                contrib[x[9]] = x[5]
        total = sum(contrib.values())
        src = 'merged' if any(x[6] == 'merged' for x in cl) else 'unmerged'
        etype = 'primary' if prim else 'secondary'
        # event_id is piped through unchanged (assigned once upstream)
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(
            rep[0], rep[1], rep[2], rep[3], rep[4], total, src,
            rep[7], rep[8], rep[9], etype))


if __name__ == '__main__':
    main()
