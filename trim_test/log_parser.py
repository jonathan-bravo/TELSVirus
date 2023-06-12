#!/usr/bin/env python3

from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--logfile')
    parser.add_argument('--outfile')
    return parser.parse_args()

def main():
    args = parse_args()
    with open(args.outfile, 'w') as o:
        with open(args.logfile, 'r') as log:
            for line in log:
                if 'sequence too short' in line:
                    next(log)
                    read_id = next(log).strip()
                    seq = next(log).strip()
                    if '_' in read_id:
                        o.write(f'{read_id}\n')
                        o.write(f'{seq}\n')
    #[print(l) for l in open(args.logfile)]

if __name__ == '__main__':
    main()