# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

| Config Variable | Default Value |
| - | - |
| `email` | "" |
| `run_id` | "test" |
| `host_organism` | "Sus Scrofa" |
| `host_genome` | "" |
| `viral_genomes` |  "resources/test/allvirusgenomes_poly-a_removed_hand_edit.fasta" |
| `strain_db` | "resources/test/strain_db.tsv" |
| `barcodes` | "resources/test/RapidBarcode.fasta" |
| `reads` | "resources/test/reads" |
| `similarity_threshold` | 0.9 |
| `crop_len` | 37 |
| `sftclp_cutoff` | 0.5 |
| `rvhaplo_max_alignments` | 75000 |
| `rvhaplo_max_alignment_bases` | 1000000000 |
| `blat_fast_map` | true |

> *NOTE: All files can exist outside of the TELSVirus directory as long as paths are correct. Alternatively files can be symbolically linked or copied to the desired location.*

## Entrez Email

TELSVirus makes requests to NCBI host reference files and for serotype, strain, and segment information so an email must be specified. If not specified the workflow will throw an error.

## Run ID

This value is used to segregate different sample runs, and could be the same for all runs as long as all samples have a different name. Howver, if different settings are used, or multiple sample runs are done, it is wise to keep results seperated.

## Host Reference Genome 

If no location for a `host_genome` is provided, TELSVirus will use the `host_organism` value to download one from NCBI. 

## Viral Reference Genomes

It is important that the `viral_genomes` fasta file be structured as follows:

```
>HV235472.1 |JP 2009213495-A/7: Recombinant Porcine Adenovirus Vector
CATCATCAATAATATACCGCACACTTTTATTGCCCCTTTTGTGGCGTGGTGATTGGCGGA
GAGGGTTGGGGGCGGCGGGCGGTGATTGGTGGAGAGGGGTGTGACGTAGCGTGGGAACGT
...
>FW304282.1 |METHODS AND COMPOSITIONS FOR INCREASING TISSUE TROPISM OF RECOMBINANT ADENOVIRAL VECTORS
TATAAACCAGTTCCACCATGGGACCGAAGAAGCAGAAGCGCGAGCTCCCCGAGGACTTCG
ATCCAGTCTACCCCTATGACGCCCCGCAGCTGCAGATCAATCCACCCTTCGTCAGCGGGG
...
>BD080521.1 |Recombinant porcine adenovirus vector
CATCATCAATAATATACCGCACACTTTTATTGCCCCTTTTGTGGCGTGGTGATTGGCGGA
GAGGGTTGGGGGCGGCGGGCGGTGATTGGTGGAGAGGGGTGTGACGTAGCGTGGGAACGT
...
```

The TELSVirus workflow uses the accession numbers at the beginning of the FASTA
read id for NCBI queries. If no location for a `strain_db` is provided, TELSVirus
will use the accessions to generate one for you. An example of this can be found in `resources/test/allvirusgenomes_poly-a_removed_hand_edit.fasta`.

## Barcodes File

This is a file of barcodes used in Nanopore sequenceing and is used for read trimming and chimera detection. An example of this file can be found in `resources/test/RapidBarcodes.fasta`.

## Crop Length

The default number of bases to be trimmed from the start and stop of an input read, value will change based on bait location in the read.

## Input Data and Test Samples

The input data specified at `reads` should be a directory that containes many sample directories - each containing one to many FASTQ files. An example of this can be seen here: `resources/test/reads/` with the provided positive and negative controls.

## Similarity Threshold

This threshold controls length clustering and duplicate filtering. A BLAT alignment qualifies when `matches + repMatches` covers at least this fraction of both the query and target lengths. Self-alignments are excluded. Each qualifying pair is ordered lexicographically by read ID, with the larger ID marked as a duplicate candidate; the existing pruning step resolves overlapping pairs. This is not random selection or a requirement for two directional alignments.

Older configs must rename `silimarity_threshold` to `similarity_threshold`.

## BLAT Fast Map

`blat_fast_map: true` enables BLAT's `-fastMap` option (the workflow default). Set it to `false` for the more sensitive, potentially much slower search. Keep this setting fixed across a similarity sweep.

## Runtime and Logs

The local profile defaults to 16 cores. Read binning buffers sequences in memory to reduce repeated gzip writes; larger inputs therefore require more RAM. Duplicate filtering processes independent PSL files in parallel, and worker/BLAT failures stop the workflow.

After successful completion, the entire run-specific `snakemake_logs` directory is deleted, including nonempty logs. Capture console output separately if needed. Benchmark files in `snakemake_benchmarks` are retained.

## Soft Clip Cutoff

This value is the percentage threshold of a read that must be soft-clipped in an alignment for that particular alignment to be removed from further processing.

## RVHaplo Resource Cutoffs

For each sample–reference pair, RVHaplo runs only when **both** conditions hold:

```text
A <= rvhaplo_max_alignments
A * L <= rvhaplo_max_alignment_bases
```

`A` is the number of mapped primary alignments (excluding unmapped, secondary,
and supplementary records), using RVHaplo's MAPQ threshold of 0. `L` is the
reference length in bases. The product is a workload proxy, not observed aligned
bases, memory bytes, or a runtime prediction. The explicit alignment limit also
protects against expensive graph construction and clustering.

Defaults are **75,000 alignments** and **1,000,000,000 alignment-reference bases**.
Equality is allowed. Inputs with no mapped primary alignments are also skipped.
Limits must be positive integers. Raise or lower them according to available
resources and acceptable waiting time; raising a limit does not guarantee success.

These provisional defaults are informed by observed TELSVirus runs, with modest
untested headroom rather than a systematically measured maximum:

- MacBook M1 Pro: 7,079 alignments against 15,072 bp completed the RVHaplo rule
  in about 4 minutes 49 seconds with 6 allocated cores.
- Linux Desktop: four references with 4,640–65,091 alignments and lengths of
  1,060–11,477 bp completed in 56 minutes 10 seconds total with 30 allocated
  threads, GPU-assisted Medaka, and reported maximum RSS of about 4.36 GiB.
- A 303,174-alignment, 30,632-bp target spent over nine hours in conditional
  SNV filtering without completing that stage. The defaults exclude it.

The successful desktop runs used two subgraphs above 50,000 alignments and one
for smaller targets. This wrapper retains that policy (`floor(A / 25000)` above
50,000), resetting it for every reference. Subgraphs reduce downstream graph
work but do not reduce the earlier conditional SNV filtering workload. They are
an approximation, so partitioning can affect reconstructed haplotypes. These
observations do not establish a universal safe limit across hardware or datasets.

A persistent `{sample}_rvhaplo_out/rvhaplo_preflight.tsv` records each target's
counts, length, product, configured limits, decision, reason, partition count,
and allocated threads. Its `run` decision means admitted by preflight, not
verified successful reconstruction. Skipped target directories contain
`rvhaplo-skipped.flag` and are excluded from the haplotype TSV even if older
FASTA output remains. A resource skip is **not a biological negative result**.

Use a new run ID when comparing cutoff policies. Direct wrapper calls now take
seven positional arguments: `OUTDIR REF_DIR SAM_DIR SAMPLE THREADS MAX_ALIGNMENTS
MAX_ALIGNMENT_BASES`.