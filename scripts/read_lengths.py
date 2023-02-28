#!/usr/bin/env python3

import argparse
import gzip
from json import dump

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()

def get_read_lengths(reads):
    read_lengths = {}
    for read in reads:
        read_id = read.split(' ')[0][1:]
        read_len = len(next(reads))
        next(reads)
        next(reads)
        read_lengths.setdefault(read_len, []).append(read_id)
    return read_lengths

def write_out(outfile, read_lengths):
    with open(outfile, 'w') as o:
        dump(read_lengths, o)

def main():
    args = parse_args()
    reads = (line.strip() for line in gzip.open(args.infile, 'rt'))
    read_lengths = get_read_lengths(reads)
    write_out(args.outfile, read_lengths)

if __name__ == '__main__':
    main()