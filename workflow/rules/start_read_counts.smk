rule start_read_count:
    input:
        f"{OUTDIR}/{{sample}}.concat.fastq.gz"
    output:
        f"{OUTDIR}/{{sample}}.start.read.count.txt"
    conda:
        "../envs/default.yaml"
    log:
        f"{LOGDIR}/{{sample}}_read_count_snakemake.log"
    shell:
        "bash workflow/scripts/start_counts.sh {input} {output}"