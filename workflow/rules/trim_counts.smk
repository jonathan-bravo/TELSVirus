rule hard_trim_count:
    input:
        f"{OUTDIR}/{{sample}}.trimmed.log"
    output:
        f"{OUTDIR}/{{sample}}.hard.trim.count.txt"
    conda:
        "../envs/default.yaml"
    log:
        f"{LOGDIR}/{{sample}}_hard_trim_count_snakemake.log"
    shell:
        "cat {input} | "
        "awk '{{ if ($0 ~ /short/) print }}' | "
        "wc -l > {output}"


rule chimeric_count:
    input:
        f"{OUTDIR}/{{sample}}.trimmed.fastq.gz"
    output:
        f"{OUTDIR}/{{sample}}.chimeric.count.txt"
    conda:
        "../envs/default.yaml"
    log:
        f"{LOGDIR}/{{sample}}_chimer_count_snakemake.log"
    shell:
        "gzcat {input} | "
        "awk '{{ if ($0 ~ /_A/ || $0 ~ /_B/) print $0 }}' | "
        "wc -l > {output}"