#!/bin/bash
#SBATCH --job-name=TLS-virl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=24:00:00
#SBATCH --output=TELSVirus_%j.log
#SBATCH --error=TELSVirus_error_%j.log

module load snakemake
module load conda

snakemake --profiles profiles/slurm