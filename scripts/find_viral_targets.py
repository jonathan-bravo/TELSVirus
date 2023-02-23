#!/usr/bin/env python3

import argparse

## WHAT THE SUBSTRAIN THAT HAS HIGHEST COVERAGE
## INCLUDE DEPTH OF COVERAGE TO BREAK TIES

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mpileup', required = True)
    parser.add_argument('--bed', required = False)
    parser.add_argument('--outfile', required = False)
    return parser.parse_args()

def get_counts(pileup):
    covered_counts = {}
    for row in pileup:
        if int(row[1]) > 0:
            try: covered_counts[row[0]] += 1
            except KeyError: covered_counts[row[0]] = 1
    return covered_counts

def write_viral_targets(total_counts, covered_counts, outfile):
    with open(outfile, "w") as o:
        for key in covered_counts.keys():
            percent = (covered_counts[key]/int(total_counts[key])) * 100
            if percent > 80.0:
                o.write(f'{key}\t0\t{total_counts[key]}')

def main():
    args = parse_args()
    pileup = (row.split('\t') for row in open(args.mpileup))
    bed = (row.split('\t') for row in open(args.bed))
    total_counts = dict([(row[0], row[2]) for row in bed])
    covered_counts = get_counts(pileup)
    write_viral_targets(total_counts, covered_counts, args.outfile)

if __name__ == '__main__':
    main()