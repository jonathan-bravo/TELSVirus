#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', help='Directory of BLAT alignment files (.psl)', required=True)
    parser.add_argument('-o', help='Output directory', required=True)
    parser.add_argument('-s', help='Similarity threshold', type=float, required=True)
    parser.add_argument('-t', help='Number of worker processes', type=int, default=1)
    return parser.parse_args()


def read_psl(psl):
    return (q.strip().split('\t') for q in open(psl))


def skip_header(qresults):
    try:
        for _ in range(5):
            next(qresults)
    except StopIteration: pass


def high_match(q, threshold):
    total_matches = int(q[0]) + int(q[2])
    return (total_matches >= threshold * int(q[14]) # amount covered is 90% of target?
            and total_matches >= threshold * int(q[10]) # amount covered is 90% of query?
            and q[9] != q[13]) # doesn't equal self


def find_dupes(qresults, threshold):
    dupes = defaultdict(list)
    for q in qresults:
        if high_match(q, threshold):
            keeper, dup = sorted([q[9], q[13]])
            dupes[keeper].append(dup)
    return dupes


def try_remove(dupes, x):
    try: dupes.pop(x)
    except KeyError: pass


def look_through_dupes(dupes, k):
    try:
        for x in dupes[k]:
            try_remove(dupes, x)
    except KeyError: pass


def remove_doubles(dupes):
    for k in list(dupes):
        look_through_dupes(dupes, k)
    return {x for s in dupes for x in dupes[s]} # results


def write_out(outfile, dup_set):
    with open(outfile, 'w') as o:
        for dup in sorted(dup_set):
            o.write(f'{dup}\n')


def process_psl(psl_path, outdir, threshold):
    cluster = os.path.splitext(os.path.basename(psl_path))[0]
    outfile = os.path.join(outdir, f'{cluster}_dupes.txt')
    qresults = read_psl(psl_path)
    skip_header(qresults)
    dupes = find_dupes(qresults, threshold)
    cleaned_dupes = remove_doubles(dupes)
    write_out(outfile, cleaned_dupes)


def main():
    args = parse_args()
    os.makedirs(args.o, exist_ok=True)
    psl_files = [os.path.join(args.p, f) for f in os.listdir(args.p) if f.endswith('.psl')]
    worker = partial(process_psl, outdir=args.o, threshold=args.s)
    with ProcessPoolExecutor(max_workers=args.t) as executor:
        # Consume results so worker exceptions fail the workflow.
        for _ in executor.map(worker, psl_files):
            pass


if __name__ == "__main__":
    main()