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

snakemake -s TELSVirus_p1 --cluster "sbatch -A {cluster.account} -q {cluster.qos} \
-c {cluster.cpus-per-task} -N {cluster.Nodes} -t {cluster.runtime} \
--mem {cluster.mem} -J {cluster.jobname} --mail-type={cluster.mail_type} \
--mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" \
--cluster-config cluster.json --jobs 100 --latency-wait 20 \
--rerun-incomplete --use-conda

snakemake -s TELSVirus_p2 --cluster "sbatch -A {cluster.account} -q {cluster.qos} \
-c {cluster.cpus-per-task} -N {cluster.Nodes} -t {cluster.runtime} \
--mem {cluster.mem} -J {cluster.jobname} --mail-type={cluster.mail_type} \
--mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" \
--cluster-config cluster.json --jobs 100 --latency-wait 20 \
--rerun-incomplete --use-conda