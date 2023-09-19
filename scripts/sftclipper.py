#!/usr/bin/env python3

from argparse import ArgumentParser
from itertools import groupby
import pysam

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--bam')
    parser.add_argument('--outfile')
    return parser.parse_args()


def parse_cigar(cigar):
    cigar_string = []
    cig_iter = groupby(cigar, lambda c: c.isdigit())
    for g, n in cig_iter:
        cigar_string.append((int("".join(n)), "".join(next(cig_iter)[1])))
    return cigar_string


def clip_check(read, cutoff):
    cigar = parse_cigar(read.cigarstring)
    soft_clip = 0
    mapped = 0
    total = 0
    for c in cigar:
        total += c[0]
        if c[1] == 'S': soft_clip += c[0]
        if c[1] == 'M' or c[1] == 'Y': mapped += c[0]
    return soft_clip/total < cutoff and mapped >= 100


def main():
    args = parse_args()
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    removed_bam = pysam.AlignmentFile(f'{args.outfile}_REMOVED.bam', "wb", template=bamfile)
    with pysam.AlignmentFile(args.outfile, "wb", template=bamfile) as outf:
        [outf.write(read) if clip_check(read, args.cutoff) else removed_bam.write(read) for read in bamfile.fetch()]
    removed_bam.close()
    bamfile.close()



if __name__ == '__main__':
    main()