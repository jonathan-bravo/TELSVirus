#!/usr/bin/env python3

import argparse
from json import dump
from Bio import Entrez, SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--email', required=True)
    return parser.parse_args()

def get_accessions(infile):
    return [line.strip().split(' ')[0][1:] for line in open(infile) if '>' in line]

def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def gen_chunk_list(accessions):
    return list(divide_chunks(accessions, 500))

def get_sources(accessions):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', id=','.join(accessions))
    sources = [record.annotations['source'] for record in SeqIO.parse(handle, 'genbank')]
    handle.close()
    return sources

def process_chunks(accession_chunk_list):
    results = []
    [results.extend(el) for el in [get_sources(chunk) for chunk in accession_chunk_list]]
    return results

def write_out(sources, accessions, outfile):
    strain_db = {}
    [strain_db.setdefault(s, []).append(a) for s,a in zip(*[sources, accessions])]
    with open(outfile, 'w') as f:
        dump(strain_db, f)

def main():
    args = parse_args()
    Entrez.email = args.email
    accessions = get_accessions(args.infile)
    accession_chunk_list = gen_chunk_list(accessions)
    sources = process_chunks(accession_chunk_list)
    write_out(sources, accessions, args.outfile)

if __name__ == '__main__':
    main()