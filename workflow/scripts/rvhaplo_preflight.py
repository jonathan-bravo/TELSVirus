#!/usr/bin/env python3
"""Measure RVHaplo inputs and emit one tab-separated preflight decision."""
import argparse
from pathlib import Path

from Bio import SeqIO
import pysam


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError('limits must be positive integers')
    return value


def decide(alignments, length, max_alignments, max_alignment_bases):
    product = alignments * length
    reasons = []
    if alignments == 0:
        reasons.append('no_mapped_primary_alignments')
    if alignments > max_alignments:
        reasons.append('max_alignments')
    if product > max_alignment_bases:
        reasons.append('max_alignment_bases')
    return product, 'skip' if reasons else 'run', ','.join(reasons) or 'within_limits'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sam', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--max-alignments', required=True, type=positive_int)
    parser.add_argument('--max-alignment-bases', required=True, type=positive_int)
    args = parser.parse_args()
    reference = SeqIO.read(args.reference, 'fasta')  # Exactly one reference per target.
    length = len(reference.seq)
    if length == 0:
        raise ValueError('Empty reference sequence')
    if not Path(args.sam).is_file():
        raise FileNotFoundError(args.sam)
    with pysam.AlignmentFile(args.sam, 'r') as sam:
        if sam.nreferences != 1 or sam.lengths[0] != length:
            raise ValueError('SAM reference length/count does not match target FASTA')
        # RVHaplo uses MAPQ 0; exclude unmapped, secondary and supplementary records.
        alignments = sum(1 for read in sam if not (read.flag & 0x904)
                         and read.mapping_quality >= 0)
    product, decision, reason = decide(
        alignments, length, args.max_alignments, args.max_alignment_bases)
    print(alignments, length, product, decision, reason, sep='\t')


if __name__ == '__main__':
    main()
