#!/usr/bin/env python3

import argparse
from collections import defaultdict
import gzip


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()


def gen_read_length_clusters(reads, outdir):
    """Group reads by sequence length, writing one gzip file per length bin."""
    buffers = defaultdict(list)  # filepath -> list of (header_bytes, seq_bytes)
    for read in reads:
        read_id = read.split(b' ')[0][1:].strip()
        seq = next(reads)
        seq_len = len(seq.strip())
        next(reads)  # desc
        next(reads)  # qual
        outfile = f'{outdir}/{seq_len}_rl_bins.fasta.gz'
        buffers[outfile].append((b'>' + read_id + b'\n', seq))

    for filepath, entries in buffers.items():
        with gzip.open(filepath, 'wb') as o:
            for header, seq in entries:
                o.write(header)
                o.write(seq)


def main():
    args = parse_args()
    reads = (line for line in gzip.open(args.infile, 'rb'))
    gen_read_length_clusters(reads, args.outdir)


if __name__ == "__main__":
    main()