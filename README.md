# README

## Cleaning the data

- [ ] Make sure all required outputs are in the `all_rule`
- [ ] Convert RVHaplo viral accession to Viral Name (+ accession ?)
- [ ] Investiagte Strainline errors
- [ ] Investigate RVHaplo errors (exit error 0)
- [ ] Check trimmed fastq data for **exact match** adapters
  - Make sure the adapter pairs are together
  - Make sure the same barcode is present and same for read from specific sample
- [ ] Check for chimeric reads
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600009/ -> Chimeric reads seem to appear in 1.7% of sequence even when samples are prepped separately
  - _if_ read has exact match of adapter align to viral reference
  - _if_ read has $\geq$ 80% match **not** chimeric
  - _else_ read **is** chimeric and is split at adapter
  - Should happen after deduplication/ host removal
    - to align we need to have already selected the target viral strain (variant)
    - if a read is chimeric then the two resulting smaller reads may belong to different viral variants that selected at start though?
    - could also affect deduplication because technically they are two different reads?
    - could also affect trimming, because if read is chimeric that could mean that the adapter sequence was not at the beginning of the sequence like we might expect?
    - So does it makes more sense to remove chimeras during read trimming? Instead of comparing to a selected viral ref, we would just compare to the whole viral database?
      - Should happen after deduplication/ host removal; but to align we need to have already selected the target viral strain (variant), but if a read is chimeric then the two resulting smaller reads may belong somewhere else? This could also affect deduplication because technically they are two different reads. Could also affect trimming, becasue if read is chimeric that could mean that the adapter sequence was not at the begging of the sequence like we might expect?

## Generate improved stats

- [ ] Write r-script(s) for looking at read qc data

## Extra

- Command to make `dag.pdf`
  - `snakemake --forceall --rulegraph | dot -Tsvg > dag.svg`

**Putting haplotype generation on hold FOR NOW**

- [ ] Add check for number of reads aligned to single strain and if (for some reason) the read count falls below a selected threshold (look at RVHaplo docs to confirm) don't run RVHaplo for this strain X barcode

## Workflow Image

![Workflow Image](dag.svg)