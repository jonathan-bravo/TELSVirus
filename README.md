# README

## Cleaning the data

- [ ] Create custom adapter file for interrogating the trimmed sequences
- [ ] Make sure that data being trimmed is actually adapter sequence
- [ ] * MAYBE change the clustering from `+40` (+- 20 bp) to `+(ceil(len(sequence)*0.11)*2)` (+- 11% bp)
  - So this will actually increase the "length" of bins that will be included for larger read lengths
  - Hopefully this will increase the number of clusters, reducing the number of reads in each cluster, and reducing the amount of time the all_vs_all alignment takes

## Generate improved stats

- [ ] Write r-script(s) for looking at read qc data

## Extra

- [x] Create DAG flowchart and include in README
  - `snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf`

**Putting haplotype generation on hold FOR NOW**

- [ ] Add check for number of reads aligned to single strain and if (for some reason) the read count falls below a selected threshold (look at RVHaplo docs to confirm) don't run RVHaplo for this strain X barcode

## Workflow Image

![Workflow Image](dag.svg)