#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from os import listdir
import pandas as pd


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--indir')
    parser.add_argument('--strains', required=True)
    parser.add_argument('--outfile')
    # Add virus db as input
    return parser.parse_args()


def parse_strains(strain_db):
    strains = (row.split('\t') for row in open(strain_db))
    return dict([(row[0], row[1]) for row in strains])


def get_strain_name(strain):
    try: name = strains[strain]
    except KeyError: name = '*'
    return name.strip()


def get_haplotypes(fq):
    return [seq for seq in SeqIO.parse(open(fq), 'fasta')]


def make_df():
    global df
    columns = [ # add virus name after strain
        'Strain',
        'Name',
        'Haplotype',
        'Length',
        'Abundance',
        'Reads#',
        'Depth',
        'Seq'
    ]
    df = pd.DataFrame(columns=columns)


def parse_haplotype(strain, hap): # Add name here after strain
    data = str(hap.id).split('_')
    df.loc[len(df.index)] = [
        strain,
        get_strain_name(strain),
        data[1],
        data[3],
        data[5],
        data[9],
        data[11],
        str(hap.seq)
    ]


def parse_rvhaplo_out(indir, directory):
    strain = directory.split('_')[-1]
    f = f'{indir}{directory}/rvhaplo_haplotypes.fasta'
    fq = get_haplotypes(f)
    [parse_haplotype(strain, hap) for hap in fq]


def main():
    global strains

    args = parse_args()

    # haplotype_0_length_979_abundance_1_number_of_reads_94_depth_50.65580448065173

    make_df()

    strains = parse_strains(args.strains)

    # pulling data from haplotype fasta files
    [parse_rvhaplo_out(args.indir, directory) for directory in listdir(args.indir)]

    df.to_csv(args.outfile, sep='\t', index=False)
        

if __name__ == '__main__':
    main()