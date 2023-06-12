#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mpileup', required = True)
    parser.add_argument('--bed', required = True)
    parser.add_argument('--strains', required = True)
    parser.add_argument('--logfile', required = True)
    parser.add_argument('--outfile', required = True)
    return parser.parse_args()

def read_pileup(pileup):
    return (row.split('\t') for row in open(pileup))

def parse_bed(bed_file):
    bed = (row.split('\t') for row in open(bed_file))
    return dict([(row[0], row[2]) for row in bed])

def parse_strains(strain_db):
    strains = (row.split('\t') for row in open(strain_db))
    return dict([(row[0], row[1]) for row in strains])

def parse_row(covered_counts_dict, strains, row):
    covered_counts_dict.setdefault(row[0], [strains[row[0]].strip(), 0, 0])
    covered_counts_dict[row[0]][1] += int(row[1])
    if int(row[1]) > 0: covered_counts_dict[row[0]][2] += 1

def get_counts(pileup, strains): 
    covered_counts = {}
    [parse_row(covered_counts, strains, row) for row in pileup]
    return covered_counts

def genome_fraction(key, cc, tc):
    return (cc[key][2]/int(tc[key])) * 100

def mean_depth(key, cc, tc):
    return cc[key][1]/int(tc[key])

def high_match(key, cc, tc):
    return genome_fraction(key, cc, tc) > 80.0

def get_viral_targets(cc, tc): # write logfile here?
    global log_lines
    best_matches = {}
    log_lines = [(cc[key][0], genome_fraction(key, cc, tc), mean_depth(key, cc, tc), key) for key in cc.keys()]
    [best_matches.setdefault(cc[key][0], []).append((genome_fraction(key, cc, tc), mean_depth(key, cc, tc), key)) for key in cc.keys() if high_match(key, cc, tc)]
    [best_matches[key].sort(reverse=True) for key in best_matches.keys()]
    return [best_matches[key][0] for key in best_matches.keys()]

def write_viral_targets(targets, tc, outfile):
    with open(outfile, "w") as o:
        [o.write(f'{t[2]}\t0\t{tc[t[2]]}') for t in targets]

def write_log(logfile):
    with open(logfile,'w') as o:
        o.write(f'Strain\tHorizontalCov\tMeanDepth\n')
        [o.write(f'{line[0]}\t{line[1]}\t{line[2]}\n') for line in log_lines]

def main():
    args = parse_args()
    pileup = read_pileup(args.mpileup)
    total_counts = parse_bed(args.bed) # length of the virus
    strains = parse_strains(args.strains)
    covered_counts = get_counts(pileup, strains)
    targets = get_viral_targets(covered_counts, total_counts)
    write_viral_targets(targets, total_counts, args.outfile)
    write_log(args.logfile)

if __name__ == '__main__':
    main()