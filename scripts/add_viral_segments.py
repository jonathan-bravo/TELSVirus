#!/usr/bin/env python3

import argparse
from Bio import Entrez


class GBFeatures:
    def __init__(self, quals):
        self.quals = {}

        for d in quals:
            try: self.quals[d['GBQualifier_name']] = d['GBQualifier_value']
            except KeyError: self.quals[d['GBQualifier_name']] = 'NA'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--email', required=True)
    return parser.parse_args()


def get_serotype_and_strain(d):
    try: serotype = d['serotype']
    except KeyError: serotype = 'NA'
    try: strain = d['strain']
    except KeyError: strain = 'NA'
    return (serotype, strain)


def get_viral_info(accessions):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='xml', id=','.join(accessions))
    data = Entrez.read(handle)
    features = [GBFeatures(e['GBSeq_feature-table'][0]['GBFeature_quals']) for e in data]
    viral_info = [get_serotype_and_strain(f.quals) for f in features]
    handle.close()
    return viral_info


def write_updated_log(outfile, log, viral_info):
    with open(outfile, 'w') as out:
        out.write('Accession\tStrain\tHorizontalCov\tMeanDepth\tSerotype\tStrain\n')
        for i, line in enumerate(log[1:]):
            out.write(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{viral_info[i][0]}\t{viral_info[i][1]}\n')


def main():
    args = parse_args()
    Entrez.email = args.email
    viral_log = [line.strip().split('\t') for line in open(args.infile, 'r')]
    accessions = [v[0] for v in viral_log if v[0] != 'Accession']
    viral_info = get_viral_info(accessions)
    write_updated_log(args.infile, viral_log, viral_info)


if __name__ == '__main__':
    main()