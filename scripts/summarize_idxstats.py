#!/usr/bin/env python3

import argparse

def parse_cmdline_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    return parser.parse_args()

def write_out(idxstats, output_file):
    with open(output_file, 'w') as o:
        o.write('Strain\tNumberOfInputReads\tMapped\tUnmapped\tPercentMapped\n')
        [o.write(f'{strain}\t{idxstats[strain][0]}\t{idxstats[strain][1]}\t{idxstats[strain][2]}\t{idxstats[strain][3]}\n') for strain in idxstats]

def mapping_stats(infile):
    idxstats = (row.strip().split('\t') for row in open(infile))
    stats = {}
    for line in idxstats:
        number_of_reads = int(line[2]) + int(line[3])
        l1 = [number_of_reads, int(line[2]), int(line[3])]
        try:
            l2 = stats[line[0]]
            stats[line[0]] = [sum(x) for x in zip(*[l1,l2])]
        except KeyError:
            stats[line[0]] = l1
    return stats

def main():
    args = parse_cmdline_params()
    idxstats = mapping_stats(args.infile) # dictionary
    [idxstats[strain].append((idxstats[strain][1]/idxstats[strain][0])*100) for strain in idxstats]
    #print(idxstats)
    write_out(idxstats, args.outfile)

if __name__ == "__main__":
    main()