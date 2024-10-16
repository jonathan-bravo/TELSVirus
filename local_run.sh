#!/bin/bash

threads=$1

snakemake -c ${threads} --use-conda --rerun-incomplete --latency-wait 20