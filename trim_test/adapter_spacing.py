#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
import gzip
import re


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--fastq')
    parser.add_argument('--adapters', nargs = '+')
    parser.add_argument('--outfile')
    return parser.parse_args()


def get_reads(fq):
    return [seq for seq in SeqIO.parse(gzip.open(fq, 'rt'), 'fastq')]


def get_adapters(fq):
    return [str(fasta.seq) for fasta in SeqIO.parse(open(fq), 'fasta')]


def find_barcode_matches(b, seq):
    return [(m.start(), m.end(), len(seq), 'b') for m in re.finditer(rf'{b}', seq)]


def find_adapter_matches(a, seq):
    return [(m.start(), m.end(), len(seq), 'a') for m in re.finditer(rf'{a}', seq)]


def barcode_check(seq):
    return [b for bar in barcodes for b in find_barcode_matches(bar, seq) if bar in seq]


def adapter_check(seq):
    return [a for adp in adapters for a in find_adapter_matches(adp, seq) if adp in seq]


def process_read(read):
    barcode_matches = barcode_check(str(read.seq))
    adapter_matches = adapter_check(str(read.seq))

    if len(adapter_matches) > 0 and len(barcode_matches) > 0:
        matches = []

        for b in barcode_matches:
            matches.append(b)
        for a in adapter_matches:
            matches.append(a)

        read_info.append(matches)


def write_out():
    with open(outfile, 'w') as o:
        [o.write(f'{m}\n') for m in read_info]


def main():
    global barcodes # 24-mers
    global adapters # 8-mers
    global read_info
    global outfile
    
    args = parse_args()
    reads = get_reads(args.fastq)
    barcodes = get_adapters(args.adapters[0])
    adapters = get_adapters(args.adapters[1])
    read_info = []
    outfile = args.outfile

    [process_read(read) for read in reads]

    write_out()


if __name__ == '__main__':
    main()