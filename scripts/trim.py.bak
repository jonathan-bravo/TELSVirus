#!/usr/bin/env python3

import argparse
import gzip
from datetime import date

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--crop', type=int, required=True)
    parser.add_argument('--logfile', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()

def get_reads(fq):
    return (line for line in gzip.open(fq, 'rb'))

def make_log(log, crop):
    with open(log, 'w') as o:
        o.write(f'{date.today()}\n\n')
        o.write(f'Cropping {crop} bases from head and tail of all reads\n\n')

def make_fq(fq):
    open(fq, 'wb').close()

def write_seq_log(logfile, read_tag, seq):
    with open(logfile, 'a') as o:
        o.write(f'>{read_tag}')
        o.write(f'{seq}\n')

def write_crop_log(logfile, read_tag, seq, crop):
    with open(logfile, 'a') as o:
        o.write(f'>{read_tag}_head\n')
        o.write(f'{seq[:crop]}\n')
        o.write(f'>{read_tag}_tail\n')
        o.write(f'{seq[-crop:]}\n\n')
       
def write_read(outfile, read_id, seq, desc, qual, crop):
    with gzip.open(outfile, 'a') as o:
        o.write(read_id)
        o.write(seq[crop:-crop] + b'\n')
        o.write(desc)
        o.write(qual[crop:-crop] + b'\n')

def process_read(read, reads, crop, outfile, logfile):
    read_id = read
    seq = next(reads)
    seq_len = len(seq.strip())
    desc = next(reads)
    qual = next(reads)
    read_tag = read_id.split(b' ')[0][1:].decode('utf-8')
    if seq_len >= (2*crop)+1:
        write_read(outfile, read_id, seq.strip(), desc, qual.strip(), crop)
        write_crop_log(logfile, read_tag, seq.strip().decode('utf-8'), crop)
    else:
        write_seq_log(logfile, read_tag, seq.decode('utf-8'))

def trim_reads(reads, crop, outfile, logfile):
    [process_read(read, reads, crop, outfile, logfile) for read in reads]

def main():
    args = parse_args()
    reads = get_reads(args.infile)
    make_log(args.logfile, args.crop)
    make_fq(args.outfile)
    trim_reads(reads, args.crop, args.outfile, args.logfile)

if __name__ == '__main__':
    main()