rule host_removal_stats:
    input:
        expand(f"{OUTDIR}/{{sample}}.remove.host.samtools.idxstats", sample = SAMPLES)
    output:
        temp(f"{OUTDIR}/host.removal.stats")
    conda:
        "../envs/alignment.yaml"
    log:
        f"{LOGDIR}/host_removal_stats_snakemake.log"
    shell:
        "python workflow/scripts/samtools_idxstats.py -i {input} -o {output}"