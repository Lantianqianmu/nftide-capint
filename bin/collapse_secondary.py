#!/usr/bin/env python3
"""Collapse secondary events into a per-primary table (HIVID step 6b).

Reads a breakpoint table that carries event_id / event_type columns (the
*_breakpoints_final.txt from filter_signal.py, or the merged/combined table)
and emits ONE ROW PER PRIMARY breakpoint. All secondary events sharing the
same event_id are collapsed into an appended 'secondary_events' column:

    chr4,190122547,+,0|chr5,10002,-,0|chr12,10027,-,0

Each secondary event is 'host_chr,host_pos,host_strand' (","-separated);
multiple secondary events of the same primary are "|"-separated. Primary rows
with no secondary events get an empty column. Secondary-only event_ids (no
primary row) are dropped.

Usage:
    collapse_secondary.py <breakpoints.txt>
"""
import sys


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: collapse_secondary.py <breakpoints.txt>\n")
        sys.exit(1)

    with open(sys.argv[1]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_eid = header.index("event_id")
            i_type = header.index("event_type")
            i_hchr = header.index("host_chr")
            i_hpos = header.index("host_pos")
            i_hstrand = header.index("host_strand")
        except ValueError as e:
            sys.stderr.write("missing required column: %s\n" % e)
            sys.exit(1)

        records = []        # list of (row_fields, event_type)
        secondaries = {}    # event_id -> list of secondary strings
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_type:
                continue
            etype = f[i_type]
            eid = f[i_eid]
            records.append((f, etype))
            if etype == "secondary":
                sec = "{},{},{}".format(f[i_hchr], f[i_hpos], f[i_hstrand])
                secondaries.setdefault(eid, []).append(sec)

        out = sys.stdout
        out.write("\t".join(header + ["secondary_events"]) + "\n")
        seen = set()
        for f, etype in records:
            if etype != "primary":
                continue
            eid = f[i_eid]
            if eid in seen:
                continue
            seen.add(eid)
            sec_col = "|".join(secondaries.get(eid, []))
            out.write("\t".join(f + [sec_col]) + "\n")


if __name__ == "__main__":
    main()
