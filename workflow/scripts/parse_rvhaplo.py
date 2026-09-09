#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import listdir
from os.path import exists, isdir, join
import pandas as pd


COLUMNS = [
    'Strain',
    'Name',
    'Haplotype',
    'Length',
    'Abundance',
    'Reads#',
    'Depth',
    'Seq'
]


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--indir')
    parser.add_argument('--strains', required=True)
    parser.add_argument('--outfile')
    return parser.parse_args()


def parse_strains(strain_db):
    strains = (row.split('\t') for row in open(strain_db))
    return dict([(row[0], row[1]) for row in strains])


def get_strain_name(strain, strains):
    try:
        name = strains[strain]
    except KeyError:
        name = '*'
    return name.strip()


def get_haplotypes(fasta_path):
    """Parse haplotypes from FASTA file, return empty list if file missing or invalid."""
    if not exists(fasta_path):
        return []
    try:
        return [seq for seq in SeqIO.parse(fasta_path, 'fasta')]
    except Exception:
        return []


def parse_haplotype(strain, hap, strains):
    """Parse a single haplotype record into a dictionary."""
    data = str(hap.id).split('_')
    try:
        return {
            'Strain': strain,
            'Name': get_strain_name(strain, strains),
            'Haplotype': data[1],
            'Length': data[3],
            'Abundance': data[5],
            'Reads#': data[9],
            'Depth': data[11],
            'Seq': str(hap.seq)
        }
    except IndexError:
        # Malformed haplotype ID, skip it
        return None


def parse_rvhaplo_out(indir, directory, strains):
    """Parse RVHaplo output directory, return list of haplotype records."""
    strain = directory.split('_')[-1]
    fasta_path = join(indir, directory, 'rvhaplo_haplotypes.fasta')
    haplotypes = get_haplotypes(fasta_path)

    records = []
    for hap in haplotypes:
        record = parse_haplotype(strain, hap, strains)
        if record is not None:
            records.append(record)
    return records


def main():
    args = parse_args()

    strains = parse_strains(args.strains)

    # Collect all haplotype records
    all_records = []

    # Handle case where input directory doesn't exist or is empty
    if args.indir and exists(args.indir) and isdir(args.indir):
        try:
            directories = listdir(args.indir)
            for directory in directories:
                dir_path = join(args.indir, directory)
                if isdir(dir_path):
                    records = parse_rvhaplo_out(args.indir, directory, strains)
                    all_records.extend(records)
        except Exception as e:
            print(f"Warning: Error reading input directory: {e}")

    # Create DataFrame (will be empty if no records found)
    df = pd.DataFrame(all_records, columns=COLUMNS)

    # Always write output file (even if empty with just headers)
    df.to_csv(args.outfile, sep='\t', index=False)

    if len(df) == 0:
        print("No haplotypes found - output file contains headers only")
    else:
        print(f"Parsed {len(df)} haplotypes")


if __name__ == '__main__':
    main()