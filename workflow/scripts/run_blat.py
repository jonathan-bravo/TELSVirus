#!/usr/bin/env python3

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path
import subprocess


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--read_clusters', required=True)
    parser.add_argument('--fast_map', action='store_true', default=False)
    return parser.parse_args()


def get_read_clusters(directory):
    return sorted(str(p.resolve()) for p in Path(directory).glob('*_rl_clusters.fasta.gz'))


def run_blat(cluster, outdir, fast_map):
    cluster_id = Path(cluster).name.split('_rl_clusters.fasta.gz')[0]
    output = Path(outdir) / f'{cluster_id}.psl'
    cmd = ['blat']
    if fast_map:
        cmd.append('-fastMap')
    cmd.extend([cluster, cluster, str(output)])
    subprocess.run(cmd, check=True)


def thread(threads, clusters, outdir, fast_map):
    with ProcessPoolExecutor(max_workers=threads) as pool:
        # Consume results so nonzero BLAT exits and worker exceptions propagate.
        for _ in pool.map(run_blat, clusters, repeat(outdir), repeat(fast_map)):
            pass


def main():
    args = parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    clusters = get_read_clusters(args.read_clusters)
    if not clusters:
        raise ValueError(f'No read clusters found in {args.read_clusters}')
    thread(args.threads, clusters, args.outdir, args.fast_map)


if __name__ == '__main__':
    main()
