#!/usr/bin/env python3

from Bio import SeqIO
from argparse import ArgumentParser
import gzip


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--fastq')
    parser.add_argument('--targets')
    parser.add_argument('--outfile')
    return parser.parse_args()


def get_reads(fq):
    return (seq for seq in SeqIO.parse(gzip.open(fq, 'rt'), 'fastq'))


def get_targets(targets):
    return [target.strip() for target in open(targets, 'r')]

    
def write_out():
    with open(outfile, 'w') as o:
        [SeqIO.write(read, o, 'fastq') for read in reads if read.id in targets]


def main():
    global reads
    global targets
    global outfile
    
    args = parse_args()
    reads = get_reads(args.fastq)
    targets = get_targets(args.targets)
    outfile = args.outfile
    write_out()


if __name__ == '__main__':
    main()