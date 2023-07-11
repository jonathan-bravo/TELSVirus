#!/usr/bin/env python3

from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--logfile')
    return parser.parse_args()


def main():
    args = parse_args()
    lines = (line.strip() for line in open(args.logfile))
    with open('cut_read_ids.txt', 'w') as o:
        for line in lines:
            if line == 'sequence too short with crop lengths:':
                next(lines)
                o.write(f'{next(lines)[1:]}\n')
        #ids = [l.next().next() for line in l if line == 'sequence too short with crop lengths:']
    #print(ids)


if __name__ == '__main__':
    main()