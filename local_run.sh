#!/bin/bash

threads=$1

snakemake -s TELSVirus_p1 -c ${threads} \
--use-conda --rerun-incomplete --latency-wait 20

snakemake -s TELSVirus_p2 -c ${threads} \
--use-conda --rerun-incomplete --latency-wait 20