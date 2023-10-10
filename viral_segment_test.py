#!/usr/bin/env python3

import argparse
from Bio import Entrez


class GBFeatures:
    def __init__(self, quals):
        self.quals = {}

        for d in quals:
            self.quals[d['GBQualifier_name']] = d['GBQualifier_value']


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile')
    parser.add_argument('--email', required=True)
    return parser.parse_args()


def get_strain(d):
    try: source = d['source']
    except KeyError: source = 'NA'
    try: strain = d['strain']
    except KeyError: strain = 'NA'
    return (source, strain)


def get_sources(accessions):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='xml', id=','.join(accessions))
    data = Entrez.read(handle)
    features = [GBFeatures(e['GBSeq_feature-table'][0]['GBFeature_quals']) for e in data]
    sources = (get_strain(f.quals) for f in features)
    handle.close()
    return sources


def write_updated_log(log, sources):
    pass


def main():
    args = parse_args()
    Entrez.email = args.email
    viral_log = [line.strip().split('\t') for line in open(args.infile, 'r')]
    accessions = [v[0] for v in viral_log if v[0] != 'Accession']
    sources = get_sources(accessions)
    write_updated_log(viral_log, sources)


if __name__ == '__main__':
    main()