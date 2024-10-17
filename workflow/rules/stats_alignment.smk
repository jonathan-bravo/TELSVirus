rule align_to_viruses_for_stats_no_secondary: # FIRST SIGN OF DUPES
    input:
        viruses = VIRUSES,
        samples = f"{OUTDIR}/{{sample}}.non.host.fastq.gz"
    output:
        temp(f"{OUTDIR}/{{sample}}.stats.viruses.sam")
    conda:
        "../envs/alignment.yaml"
    threads:
        32
    benchmark:
        f"{BENCHDIR}/{{sample}}_ato_viruses_4stats.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_ato_viruses_4stats_snakemake.log"
    shell:
        "minimap2 "
        "--secondary=no "
        "--q-occ-frac 0 "
        "-f 0 "
        "-t {threads} "
        "-o {output} "
        "-a {input.viruses} {input.samples}"