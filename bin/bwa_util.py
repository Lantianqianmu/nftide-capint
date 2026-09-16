#!/usr/bin/env python3
"""Enumerate found alternative host placements, retaining full SAM geometry."""
import os
import subprocess
import tempfile
from junctions import collect


def remap_host_all(reads, host_index, min_score, threads=2):
    if not reads or not host_index:
        return {}
    with tempfile.TemporaryDirectory(prefix='hbv_alternatives_') as tmp:
        fa, sam = os.path.join(tmp,'reads.fa'), os.path.join(tmp,'host.sam')
        with open(fa,'w') as out:
            for name, seq in reads:
                out.write('>'+name+'\n'+seq+'\n')
        with open(sam,'w') as out:
            subprocess.run(['bwa','mem','-a','-Y','-T',str(min_score),'-t',str(threads),host_index,fa], stdout=out, check=True)
        return collect(sam)
