#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import Entrez
from tqdm import tqdm

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--email', required=True)
    return parser.parse_args()

def get_accessions(infile):
    return [line.strip().split(' ')[0][1:] for line in open(infile) if '>' in line]

def gen_chunk_list(accessions, n):
    return [accessions[i:i + n] for i in range(0, len(accessions), n)]

def get_sources(accessions):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', id=','.join(accessions))
    sources = [line.strip().split('     ')[-1][1:] for line in handle if 'SOURCE' in line]
    handle.close()
    return sources

def process_chunks(accession_chunk_list):
    results = []
    [results.extend(el) for el in (get_sources(chunk) for chunk in tqdm(accession_chunk_list))]
    return results

def write_out(accessions, sources, outfile):
    with open(outfile, 'w') as o:
        [o.write(f'{a}\t{s}\n') for a,s in zip(*[accessions, sources])]

def main():
    args = parse_args()
    Entrez.email = args.email
    accessions = get_accessions(args.infile)
    accession_chunk_list = gen_chunk_list(accessions, 1000)
    sources = process_chunks(accession_chunk_list)
    write_out(accessions, sources, args.outfile)

if __name__ == '__main__':
    main()