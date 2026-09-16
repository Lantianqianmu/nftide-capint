#!/usr/bin/env python3
"""Append NNSS to the evidence table without dropping uncertainty metadata."""
import argparse
import csv
from table_io import dict_reader
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('efr',type=int)
    p.add_argument('--nnss-threshold',type=float,default=1.0,help='Deprecated; no filtering')
    args=p.parse_args()
    if args.efr<0:
        p.error('EFR must be nonnegative')
    reader=dict_reader(sys.stdin)
    writer=csv.DictWriter(sys.stdout,fieldnames=list(reader.fieldnames)+['efr','nnss'],delimiter='\t',lineterminator='\n')
    writer.writeheader()
    for row in reader:
        row['efr']=args.efr
        row['nnss']='%.3f' % (int(row['nss'])*1e6/args.efr if args.efr else 0)
        writer.writerow(row)

if __name__=='__main__':
    main()
