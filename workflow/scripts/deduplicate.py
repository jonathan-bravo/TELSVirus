#!/usr/bin/env python3

import argparse
import gzip


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--duplicates', required=True)
    parser.add_argument('--out_reads', required=True)
    parser.add_argument('--out_dupes', required=True)
    return parser.parse_args()


def get_duplicates(dupes_file):
    return {line.strip() for line in open(dupes_file)}


def get_reads(reads):
    return (line for line in gzip.open(reads, 'rb'))


def write_out(reads_out, dupes_out, is_duplicate, read_id, seq, desc, qual):
    """Write a single FASTQ record to the appropriate output stream."""
    out_stream = dupes_out if is_duplicate else reads_out
    out_stream.write(read_id)
    out_stream.write(seq)
    out_stream.write(desc)
    out_stream.write(qual)


def dup_check(read, reads, duplicates, out_reads, out_dupes):
    read_id = read.split(b' ')[0][1:].strip().decode('utf-8')
    is_dup = read_id in duplicates
    write_out(out_reads, out_dupes, is_dup, read, next(reads), next(reads), next(reads))


def gen_deduped_reads(reads, duplicates, out_reads, out_dupes):
    """Stream reads through, writing to appropriate output files."""
    with gzip.open(out_reads, 'wb') as reads_out, gzip.open(out_dupes, 'wb') as dupes_out:
        for read in reads:
            dup_check(read, reads, duplicates, reads_out, dupes_out)


def main():
    args = parse_args()
    duplicates = get_duplicates(args.duplicates)
    reads = get_reads(args.reads)
    gen_deduped_reads(reads, duplicates, args.out_reads, args.out_dupes)


if __name__ == "__main__":
    main()