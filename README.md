# TELSVirus

[![Static Badge](https://img.shields.io/badge/snakemake-%E2%89%A57.12.1-%23039475?logo=data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICB3aWR0aD0iMjMuNjQ5OTQ0bW0iCiAgIGhlaWdodD0iMjMuNjUyNzM5bW0iCiAgIHZpZXdCb3g9IjAgMCAyMy42NDk5NDQgMjMuNjUyNzM5IgogICB2ZXJzaW9uPSIxLjEiCiAgIGlkPSJzdmcxIgogICB4bWw6c3BhY2U9InByZXNlcnZlIgogICB4bWxuczp4bGluaz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94bGluayIKICAgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIgogICB4bWxuczpzdmc9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48ZGVmcwogICAgIGlkPSJkZWZzMSI+PGNvbG9yLXByb2ZpbGUKICAgICAgIG5hbWU9IkRpc3BsYXkiCiAgICAgICB4bGluazpocmVmPSJmaWxlOi8vL0xpYnJhcnkvQ29sb3JTeW5jL1Byb2ZpbGVzL0Rpc3BsYXlzLy0xMUI3M0RDNy0wMDhFLTQ0MUUtQjQ4Mi02MTk3MkVENTQ5MzAuaWNjIgogICAgICAgaWQ9ImNvbG9yLXByb2ZpbGUxIiAvPjxjb2xvci1wcm9maWxlCiAgICAgICBuYW1lPSJBQ0VTLUNHLUxpbmVhci1BY2FkZW15LUNvbG9yLUVuY29kaW5nLVN5c3RlbS1BUDEiCiAgICAgICB4bGluazpocmVmPSJmaWxlOi8vL1N5c3RlbS9MaWJyYXJ5L0NvbG9yU3luYy9Qcm9maWxlcy9BQ0VTQ0clMjBMaW5lYXIuaWNjIgogICAgICAgaWQ9ImNvbG9yLXByb2ZpbGUyIiAvPjwvZGVmcz48ZwogICAgIGlkPSJsYXllcjEiCiAgICAgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTEwOS45MDcyMywtMTMxLjQ2NjUxKSI+PHBhdGgKICAgICAgIHN0eWxlPSJmaWxsOiMwMzk0NzU7ZmlsbC1vcGFjaXR5OjEiCiAgICAgICBkPSJtIDExNS41Mzg5MSwxNTQuODg4NzEgYyAtMy41NzMwNiwtMC43NzczNiAtNi41NTcwMiwtNC45OTY3NSAtNC40NDQ2LC02LjI4NDc4IDAuNTMzMTcsLTAuMzI1MDkgMTQuNzg4OTEsLTAuNDUzNjEgMTUuOTAzODcsLTAuMTQzMzcgMi44Mzc4OCwwLjc4OTYzIDMuMDk1MzYsNS4xMzQ1IDAuMzc3NDMsNi4zNjkwMiAtMC43NzE4MiwwLjM1MDU3IC0xMC4yNzg4MiwwLjM5ODA2IC0xMS44MzY3LDAuMDU5MSB6IG0gMTUuMjcyNTksLTMuMTUyMTIgYyAtMC4xOTQyNywtNC4wNjMyNiAtMS43ODYzLC00LjkyNTgzIC05LjMwMTkxLC01LjAzOTg1IC01LjI2NzA4LC0wLjA3OTkgLTUuMjY3MDgsLTAuMDc5OSAtNS45MDQ3OSwtMC41Mjg0MSAtMi4wOTM5OSwtMS40NzI3MyAtMi4wMjU2MywtNC41MTU1MiAwLjEyOTksLTUuNzgyNTcgMC45OTk3MywtMC41ODc2NiAxMS4yODYwMiwtMC42NTA2MyAxMi44MDIxNCwtMC4wNzg0IDQuOTI1MjEsMS44NTg5OCA2LjY1ODc4LDguMzYxOTUgMy4yMzg0LDEyLjE0NzkgLTAuNzg2MTMsMC44NzAxNiAtMC44OTE1NCwwLjc5MTU2IC0wLjk2Mzc0LC0wLjcxODY5IHogbSAtMTkuMzA3NzcsLTcuODkzMzYgYyAtMy4xNjA2OSwtNC4xODYyOSAtMS40MzQ2NywtMTAuMTI3NjQgMy40NzA5MSwtMTEuOTQ3NjcgMS4xODg3MSwtMC40NDEwMiAxMC43MzczMywtMC41OTY0NCAxMi42MTQyNCwtMC4yMDUzMSAzLjMwNDgyLDAuNjg4NjggNi4yMzEyNiwzLjkzMjQ2IDUuMjUzNTIsNS44MjMxOSAtMC4zODkyNywwLjc1Mjc2IC0wLjQ2MzY2LDAuNzU5NzkgLTguMDM5MTYsMC43NTk3OSAtOC43MzM1MSwwIC05LjAyNjE2LDAuMDQzNCAtMTAuNjcxOTEsMS41ODI4MSAtMS4wMzQ1OSwwLjk2Nzc0IC0xLjQ5NzI0LDIuMDAzODcgLTEuNTc3NSwzLjUzMjg5IC0wLjA3ODMsMS40OTE4MSAtMC4yMjAzNywxLjU1MzI3IC0xLjA1MDEsMC40NTQzIHogbSAxNi41MDM2MiwtOC43NTkxOCBjIDAuNDY5ODgsLTAuMzI5MTIgMC4zOTIxNywtMS4wNzYzOSAtMC4xNDMwMiwtMS4zNzUzMyAtMC44NjIwNCwtMC40ODE1MSAtMS42NDM0NiwwLjgxMzY3IC0wLjgzMzI5LDEuMzgxMTQgMC40MTIzLDAuMjg4NzggMC41NTY5NCwwLjI4NzkyIDAuOTc2MzEsLTAuMDA2IHoiCiAgICAgICBpZD0icGF0aDEiIC8+PC9nPjwvc3ZnPgo=)](https://snakemake.readthedocs.io/en/stable/)
[![Static Badge](https://img.shields.io/badge/conda-%E2%89%A54.12.0-%2344A833?logo=anaconda&logoColor=%2344A833)](https://anaconda.org/)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)
![Lint Code Base](https://github.com/jonathan-bravo/TELSVirus/actions/workflows/linter.yaml/badge.svg?event=push)

A snakemake workflow for viral strain detection.

## Requirements

All dependancies are managed through
[`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
environments included in the repository.

The workflow is currently built around basecalled Nanopore sequencing output.
This does not mean it cannot work for PacBio sequencing data, but it has not
been tested on this.

### Install Snakemake and Clone the Repository

**Create the environment for `telsvirus` using conda:**

```bash
conda create -c conda-forge -c bioconda -c anaconda -n telsvirus snakemake git git-lfs
```

**Clone the repository:**

```bash
conda activate telsvirus

git clone https://github.com/jonathan-bravo/TELSVirus.git
```

### Update Config

Instructions on updating the configuration can be found [here](config/README.md).

### Usage on Local Desktop or Interactive HPC Run

Make sure to update the `core` value `local` or `hpc` profiles located at
`workflows/profiles/local/config.yaml` or `workflows/profiles/hpc/config.yaml`
if a different number of CPU cores is available on your system.

| Profile | Profile Variable | Default Value |
| - | - | - |
| `local` | `cores` | 6 |
| `hpc` | `cores` | 120 |

**Running the workflow locally:**

```bash
cd TELSVirus

snakemake --profile worflow/profiles/local
```

**Running the workflow on an HPC interactively:**

```bash
cd TELSVirus

snakemake --profile worflow/profiles/hpc
```

### Usage on Slurm Cluster

Make sure to update the `email`, `account`, and `qos` values in the slurm
profile located at `worflow/profiles/slurm/config.yaml`

```yaml
default-resources:
  - mem_mb=32000
  - account=
  - qos=
  - email=
  - mail_type="NONE"
```

Make sure all string values are surrounded by double quotes **("")**.

**Move the `run.sh` from the `resources` directory up one level:**

```bash
mv resources/run.sh .
```

Make sure to edit the `email` and `time` if necessery for your run. *(I believe
the email is necessary for batch runs.)*

```sh
#SBATCH --mail-user=<email>
#SBATCH --time=24:00:00
```

**Launching a SLURM job for the workflow:**

```bash
cd TELSVirus

# Run the workflow
sbatch run.sh
```

## Test Data

A negative and positive sample are included in `resources/test/reads/`.

> *NOTE: `git-lfs` is a requirement for the test data to work. Without it, the*
> *FASTA and FASTQ files come through as `git-lfs` parts and will cause the*
> *workflow to error out.*

## Output

| Name | Content |
| - | - |
| `on_target_stats.tsv` | A file that containes a row with the number of input reads, number of reads mapped to host, number of reads mapped to viral database, number of unmapped reads, the host reads percent, and the on-target percent	for each sample. |
| `{sample}_add_sample_info.done` | A flag file ensureing run id and sample id are added to all metric files. |
| `{sample}_chimeric_count.txt` | A file that contains a single count of reads that were split as chimeras during trimming. |
| `{sample}_dedup.fastq.gz` | The deduplicated input reads. |
| `{sample}_dup_reads.fastq.gz` | The reads removed during deduplication. |
| `{sample}_duplicates.txt` | The ids of reads considered duplicates. |
| `{sample}_find_duplcates.done` | A flag file ensuring deduplication is finished. |
| `{sample}_hard_trim_count.txt` | Reads that were removed from analysis for being too short. |
| `{sample}_non_host.fastq.gz` | Deduplicated and host removed reads. |
| `{sample}_post_dedup_rl.tsv` | A file that contains read lengths before deduplication. |
| `{sample}_pre_dedup_rl.tsv` | A file that contains read lengths after deduplication. |
| `{sample}_reads_per_strain_filtered.tsv` | The number of reads that aligned to each viral strain in the `viral_genomes` filtered to only those with $>0$ reads. |
| `{sample}_reads_per_strain.tsv` | The number of reads that aligned to each viral strain in the `viral_genomes`. |
| `{sample}_selected_viral_targets.log` | A file that contains the selected viral strains from `viral_genomes`. A strain is selected if it has a horizontal coverage of $\ge 80\%$. If there are multiple viral accessions with the same strain then the highest horizontal coverage is chosen. If the horizontal coverage is the same then the accession with the highest mean depth is chosen.  |
| `{sample}_start_read_count.txt` | A file that contains a single count of reads before any processing. |
| `{sample}_stats_viruses_sorted_sftclp_REMOVED.bam` | Alignments that were removed from the `{sample}_stats_viruses_sorted_sftclp.bam` for failing the soft-clip check. |
| `{sample}_trimmed.fastq.gz` | The reads after trimming. |
| `{sample}_trimmed.log` | A log file containing the sequences trimmed from each read, the number of bases trimmed from each end, and the full sequence if it was too short and removed from further processing. |
| `{sample}_viral_target_genomes.fasta` | A FASTA file contining all viral sequences that are found in the `{sample}_selected_viral_targets.log` file. |
| `{sample}_VIRAL_TARGETS_FOUND` **OR** `{sample}_NO_VIRAL_TARGETS` | A flag file indicating in viral targets are found. Originally used for further processing; currently just for information. |
| `{sample}_viral_targets.log` | All viral targets before applying filtering. |
| `{sample}_viruses_sorted_sftclp_REMOVED.bam` | Alignments that were removed from the `{sample}_viruses_sorted_sftclp.bam` for failing the soft-clip check. |
| `{sample}_viruses_sorted_sftclp.bam` | The alignment files used for determing all viral targets. |
| `{sample}_viruses_sorted_sftclp.bam.bai` | The index file of `{sample}_viruses_sorted_sftclp.bam`. |
| `{sample}.mpileup` | Pileup file generated for determining viral targets horizontal coverage and mean depth. |

## Making a Workflow DAG

```bash
snakemake --forceall --rulegraph | dot -Tsvg > dag.svg
```

### Workflow DAG Image

![Workflow Image](resources/dag.svg)