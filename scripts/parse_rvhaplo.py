#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from os import listdir
import pandas as pd


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--indir')
    parser.add_argument('--outfile')
    return parser.parse_args()


def get_haplotypes(fq):
    return [seq for seq in SeqIO.parse(open(fq), 'fasta')]


def make_df():
    global df
    columns = [
        'Strain',
        'Haplotype',
        'Length',
        'Abundance',
        'Reads#',
        'Depth',
        'Seq'
    ]
    df = pd.DataFrame(columns=columns)


def parse_haplotype(strain, hap):
    data = str(hap.id).split('_')
    df.loc[len(df.index)] = [
        strain,
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
    args = parse_args()

    # haplotype_0_length_979_abundance_1_number_of_reads_94_depth_50.65580448065173

    make_df()

    # pulling data from haplotype fasta files
    [parse_rvhaplo_out(args.indir, directory) for directory in listdir(args.indir)]

    df.to_csv(args.outfile, sep='\t', index=False)
        

if __name__ == '__main__':
    main()