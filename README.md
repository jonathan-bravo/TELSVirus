# README

## Cleaning the data

- [ ] Add read deduplication
- [ ] Add read trimming using `trimmomatic`
- [ ] Create custom adapter file for read trimming

## Generate improved stats

- [x] Update alignment to host (swine) to include `--secondary=no`
  - This will improve the alignments stats removing the ~9% read inflation
- [x] Add additional alignment step with `--secondary=no` to viral db for read stats
- [x] Create an on-target reads summary -> Total Reads, Reads Mapped to Host, Reads Mapped to All Viral DB
- [x] Use alignment WITH `--secondary=yes` to generate read stats for all strains per sample
- [ ] Write script to summarize per strain idxstats

## Improving target strain selection

- [ ] Create 'Strain DB' to associate fasta header with correct virus
- [ ] Update `find_viral_targets.py` to select best strain for each virus with genome fraction >= 80% (using average depth to break ties)
  - Make special case for influenza as the virus has multiple segments for each strain

## Extra

**Putting haplotype generation on hold FOR NOW**

- [ ] Reach out to HiPerGator about symlink with daccord and Strainline
- [x] Update the alignment to single viral strains to inlcude `--secondary=no`
- [ ] Add check for number of reads aligned to single strain and if (for some reason) the read count falls below a selected threshold (look at RVHaplo docs to confirm) don't run RVHaplo for this strain X barcode