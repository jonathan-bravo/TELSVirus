#!/usr/bin/env python3

from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--sample_id')
    parser.add_argument('--run_id')
    parser.add_argument('--infiles', nargs='+')
    return parser.parse_args()


def main():
    sample_id = 'barcode01'
    run_id = 'June_07_2023_SHIC_Plate1_3'
    args = parse_args()
    for f in args.infiles:
        original = [line for line in open(f,'r')]
        with open(f, 'w') as final:
            final.write(f'run_id\tsample_id\t{original[0]}')
            for l in original[1:]:
                final.write(f'{run_id}\t{sample_id}\t{l}')


if __name__ == '__main__':
    main()