#!/usr/bin/env python3
"""BBMerge wrapper implementing the Nextflow PE_ASSEMBLE output contract.

Outputs <prefix>_merged.fq.gz, <prefix>_unmerged_R1.fq.gz and
<prefix>_unmerged_R2.fq.gz plus assembly log/QC. BBMerge supplies all merged
sequences and qualities. Only mate suffixes are normalized: none for merged
reads, /1 and /2 for unmerged mates. Molecule/copy tags remain intact.
"""
import argparse
from collections import Counter
import gzip
import itertools
import json
from pathlib import Path
import subprocess
from filter_chimeric import fastq


def paired_records(first, second):
    for a,b in itertools.zip_longest(first,second):
        if a is None or b is None or a[0]!=b[0]:
            raise ValueError('FASTQ mates are missing or out of order')
        yield a,b


def write_record(out,record,suffix=''):
    name,seq,qual=record
    out.write('@'+name+suffix+'\n'+seq+'\n+\n'+qual+'\n')


def validate_counts(expected,merged,unmerged):
    observed=merged+unmerged
    if observed!=expected:
        missing=sum((expected-observed).values())
        extra=sum((observed-expected).values())
        raise ValueError(f'BBMerge did not preserve input pairs: missing={missing}, extra={extra}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('r1');p.add_argument('r2');p.add_argument('--prefix',required=True)
    p.add_argument('--min-overlap',type=int,default=6)
    p.add_argument('--max-mismatch',type=float,default=.2)
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--keep-unmerged',choices=('true','false'),default='true')
    args=p.parse_args()
    if args.threads < 1 or args.min_overlap < 1 or not 0 <= args.max_mismatch <= 1:
        p.error('threads and min-overlap must be positive; max-mismatch must be between 0 and 1')
    prefix=args.prefix
    interleaved=prefix+'_input_interleaved.fq.gz'
    merged=prefix+'_bbmerge_merged.fq.gz'
    unmerged=prefix+'_unmerged.fq.gz'
    log=prefix+'_assembly.log'
    expected=Counter()
    with gzip.open(interleaved,'wt') as out:
        for a,b in paired_records(fastq(args.r1),fastq(args.r2)):
            expected[a[0]]+=1
            write_record(out,a,'/1');write_record(out,b,'/2')
    cmd=['bbmerge.sh','in='+interleaved,'interleaved=t','out='+merged,'outu='+unmerged,
         'minoverlap='+str(args.min_overlap),'minoverlap0='+str(args.min_overlap),
         'maxratio='+str(args.max_mismatch),'threads='+str(args.threads)]
    with open(log,'w') as out:
        subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,check=True)
    messages=Path(log).read_text()
    if 'Exception in thread' in messages or 'AssertionError' in messages:
        raise RuntimeError('BBMerge reported a worker failure; see '+log)
    merged_names=Counter(r[0] for r in fastq(merged))
    unmerged_names=Counter()
    reads=iter(fastq(unmerged))
    for a,b in paired_records(reads,reads):
        unmerged_names[a[0]]+=1
    validate_counts(expected,merged_names,unmerged_names)

    # Validate native outputs first; stage final files before exposing them.
    # Do not reconstruct consensus sequences or replace BBMerge qualities.
    outputs=[prefix+'_merged.fq.gz',prefix+'_unmerged_R1.fq.gz',prefix+'_unmerged_R2.fq.gz']
    temporary=[Path(name+'.tmp') for name in outputs]
    try:
        with gzip.open(temporary[0],'wt') as out:
            for record in fastq(merged):
                write_record(out,record)
        with gzip.open(temporary[1],'wt') as out1,gzip.open(temporary[2],'wt') as out2:
            if args.keep_unmerged=='true':
                reads=iter(fastq(unmerged))
                for a,b in paired_records(reads,reads):
                    write_record(out1,a,'/1');write_record(out2,b,'/2')
        for temporary_path,destination in zip(temporary,outputs):
            temporary_path.replace(destination)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
    summary=dict(input_pairs=sum(expected.values()),merged_pairs=sum(merged_names.values()),
                 unmerged_pairs=sum(unmerged_names.values()),keep_unmerged=args.keep_unmerged=='true',
                 emitted_unmerged_pairs=sum(unmerged_names.values()) if args.keep_unmerged=='true' else 0)
    Path(prefix+'_assembly_qc.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__=='__main__':main()
