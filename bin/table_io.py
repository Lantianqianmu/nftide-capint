"""Read quoted breakpoint TSVs, including large JSON evidence fields.

CSV's 128-KiB default is too small for repeat-rich integration candidates.
Use the largest practical platform-supported limit without truncating evidence
or changing the existing quoting format (including cached Nextflow outputs).
"""
import csv
import sys


def dict_reader(handle):
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            # Some Python/platform combinations accept a narrower C integer.
            limit //= 10
    return csv.DictReader(handle, delimiter='\t')
