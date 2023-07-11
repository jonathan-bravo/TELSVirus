#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--fq')
    parser.add_argument('--log')
    return parser.parse_args()

def main():
    args = parse_args()

if __name__ == '__main__':
    main()