#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
import gzip

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--fq')
    parser.add_argument('--outfile')
    return parser.parse_args()

def get_reads(fq):
    return (seq for seq in SeqIO.parse(gzip.open(fq, 'rt'), 'fastq') if '_' in seq.id)

def write_out(outfile, reads):
    with open(outfile, 'w') as o:
        for read in reads:
            o.write(f'>{read.id}\n')
            o.write(f'{read.seq}\n')

def main():
    args = parse_args()
    reads = get_reads(args.fq)
    write_out(args.outfile, reads)

if __name__ == '__main__':
    main()