#!/usr/bin/env python3
"""Emit one row per candidate, including ambiguous candidates and inline alternatives.

Alternative placements are already stored in secondary_events JSON and never
count as extra biological events. The historical primary filename is retained.
"""
import argparse
import csv
from table_io import dict_reader
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('breakpoints')
    args=p.parse_args()
    with open(args.breakpoints) as handle:
        reader=dict_reader(handle)
        writer=csv.DictWriter(sys.stdout,fieldnames=reader.fieldnames,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        for row in reader:
            writer.writerow(row)

if __name__=='__main__':
    main()
