#!/usr/bin/env python3

import argparse
import gzip
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--duplicates', required=True)
    parser.add_argument('--out_reads', required=True)
    parser.add_argument('--out_dupes', required=True)
    return parser.parse_args()

def get_duplicates(dupes_file):
    return [line.strip() for line in open(dupes_file)]

def get_reads(reads):
    return (line for line in gzip.open(reads, 'rb'))

def write_out(outfile, read_id, seq, desc, qual):
    with open(outfile , "a") as o:
        o.write(read_id + b'\n')
        o.write(seq + b'\n')
        o.write(desc + b'\n')
        o.write(qual + b'\n')

def gen_deduped_reads(reads, duplicates, out_reads, out_dupes):
    for read in reads:
        read_id = read.split(b' ')[0][1:]
        if read_id in duplicates:
            write_out(out_dupes, read, next(reads), next(reads), next(reads))
        else:
            write_out(out_reads, read, next(reads), next(reads), next(reads))

def main():
    args = parse_args()
    duplicates = get_duplicates(args.duplicates)
    reads = get_reads(args.reads)
    gen_deduped_reads(reads, duplicates, args.out_reads, args.out_dupes)

if __name__ == "__main__":
    main()