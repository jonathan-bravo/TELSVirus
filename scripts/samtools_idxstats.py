#!/usr/bin/env python3

import argparse
import os

def parse_cmdline_params():
    info = "Parses a Samtools idxstats file to obtain the total number of mapped, unmapped, and total reads"
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument(
        '-i', '--input_files', nargs='+', required=True,
        help='Use globstar to pass a list of files, (Ex: *.tsv)')
    parser.add_argument(
        '-o', '--output_file', required=True,
        help='Output file to write mapping results to')
    return parser.parse_args()

def header(output_file):
    with open(output_file, 'a') as o:
        o.write('Sample\tNumberOfInputReads\tMapped\tUnmapped\n')

def mapping_stats(input_list, output_file):
    for f in input_list:
        mapped = 0
        unmapped = 0
        number_of_reads = 0
        with open(f, 'r') as fp:
            sample_name = os.path.basename(str(fp.name)).split('.', 1)[0]
            for line in fp:
                columns = line.strip().split('\t')
                mapped += int(columns[2])
                unmapped += int(columns[3])
            number_of_reads += mapped + unmapped
        with open(output_file, 'a') as o:
            o.write(f'{sample_name}\t{number_of_reads}\t{mapped}\t{unmapped}\n')

def main():
    opts = parse_cmdline_params()
    header(opts.output_file)
    mapping_stats(opts.input_files, opts.output_file)

if __name__ == "__main__":
    main()