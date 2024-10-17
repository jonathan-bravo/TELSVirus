rule concat_parts:
    input:
        f"{READS}/{{sample}}"
    output:
        temp(f"{OUTDIR}/{{sample}}.concat.fastq.gz")
    conda:
        "../envs/default.yaml"
    log:
        f"{LOGDIR}/{{sample}}_concat_parts_snakemake.log"
    shell:
        "cat {input}/* > {output}"