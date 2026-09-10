#!/usr/bin/env python3

import argparse
from math import ceil
from os import listdir
from pathlib import Path
import gzip
import shutil

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', type = float, required = True)
    parser.add_argument('--indir', required = True)
    parser.add_argument('--outdir', required = True)
    return parser.parse_args()

def get_bins(indir):
    bins = [int(f.split('_')[0]) for f in listdir(indir) if f != '.DS_Store']
    bins.sort()
    return bins

def get_clusters(bins, threshold):
    val = 1.99 - threshold
    clusters = []
    while bins:
        cluster_len = ceil(bins[0]*val)
        cluster = [x for x in bins if x <= cluster_len]
        clusters.append(cluster)
        bins = [x for x in bins if x not in cluster]
    return clusters

def cat_files(clusters, indir, outdir):
    """Stream compressed bins into plain FASTA for BLAT (no gzip subprocesses)."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for cluster in clusters:
        start = cluster[0]
        end = cluster[-1]
        outfile = Path(outdir) / f'{start}_to_{end}_rl_clusters.fasta'
        with outfile.open('wb') as output:
            for length in cluster:
                infile = Path(indir) / f'{length}_rl_bins.fasta.gz'
                with gzip.open(infile, 'rb') as source:
                    shutil.copyfileobj(source, output)

def main():
    args = parse_args()
    bins = get_bins(args.indir)
    clusters = get_clusters(bins, args.threshold)
    cat_files(clusters, args.indir, args.outdir)

if __name__ == '__main__':
    main()