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

This value determines how similar two sequences have to be to eachother to be considered duplicates. If the query read (A) matches against a target read (B) at 90% and read B matches against A at 90% then one read is randomly removed from further processing.

## Soft Clip Cutoff

This value is the percentage threshold of a read that must be soft-clipped in an alignment for that particular alignment to be removed from further processing.