#!/usr/bin/env python3

import argparse
import gzip
from datetime import date
from Bio import SeqIO
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--crop', type=int)
    parser.add_argument('--logfile')
    parser.add_argument('--outfile')
    parser.add_argument('--adapters', nargs='+')
    return parser.parse_args()

def get_reads(fq):
    return (seq for seq in SeqIO.parse(gzip.open(fq, 'rt'), 'fastq'))

def get_adapters(fq):
    return [str(fasta.seq) for fasta in SeqIO.parse(open(fq), 'fasta')]

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

# def process_read(read, reads, crop, outfile, logfile):
#     read_id = read
#     seq = next(reads)
#     seq_len = len(seq.strip())
#     desc = next(reads)
#     qual = next(reads)
#     read_tag = read_id.split(b' ')[0][1:].decode('utf-8')
#     if seq_len >= (2*crop)+1:
#         write_read(outfile, read_id, seq.strip(), desc, qual.strip(), crop)
#         write_crop_log(logfile, read_tag, seq.strip().decode('utf-8'), crop)
#     else:
#         write_seq_log(logfile, read_tag, seq.decode('utf-8'))

def in_head(seq, span):
    return span[1] <= len(seq)//2

def in_tail(seq, span):
    return span[0] >= len(seq)//2

def centered(seq, span):
    return span[0] < len(seq)//2 and span[1] > len(seq)//2

def process_read(read, adapters, barcodes, crop, outfile, logfile):
        # if any(barcode in read for barcode in barcodes):
        #     print('yes')
        head_crop = crop
        tail_crop = crop
        not_found = True
        seq = str(read.seq)
        for barcode in barcodes:
            if barcode in seq:
                not_found = False
                span = re.search(rf'{barcode}', seq).span()
                all_match = re.findall(rf'{barcode}', seq)
                matches = []
                for match in re.finditer(rf'{barcode}', seq):
                    matches.append((match.start(), match.end(), match.group()))
                print(read.id, matches)
        if not_found:
            print(read.id, 'NONE')
        #         print(all_match)
        #         if centered(seq, span):
        #             print('center')
        #             #print(f'in {read.id} len {len(seq)} center {span}')
        #         else:
        #             if in_head(seq, span):
        #                 #head_crop = pos + 24
        #                 print('head')
        #                 #print(f'in {read.id} len {len(seq)} head {pos}')
        #             if in_tail(seq, span):
        #                 print('tail')
        #                 #print(f'in {read.id} len {len(seq)} tail {pos}')
        # # if not_found:
        # #     print(f'just crop {read.id}')



# def trim_reads(reads, adapters, barcode, crop, outfile, logfile):
#     [process_read(read, reads, crop, outfile, logfile) for read in reads]

def trim_reads(reads, adapters, barcodes, crop, outfile, logfile):
    [process_read(read, adapters, barcodes, crop, outfile, logfile) for read in reads]

def main():
    args = parse_args()
    reads = get_reads(args.infile)
    illumina_adapters = get_adapters(args.adapters[0])
    barcodes = get_adapters(args.adapters[1])
    #for seq in reads: print(seq)
    make_log(args.logfile, args.crop)
    make_fq(args.outfile)
    trim_reads(reads, illumina_adapters, barcodes, args.crop, args.outfile, args.logfile)

if __name__ == '__main__':
    main()