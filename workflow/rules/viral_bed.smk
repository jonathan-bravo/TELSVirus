rule viral_bed:
    input:
        VIRUSES,
    output:
        temp(f"{OUTDIR}/all.viral.targets.bed"),
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/viral_bed.benchmark"
    log:
        f"{LOGDIR}/viral_bed_snakemake.log",
    shell:
        "samtools faidx {input}; "
        "cat {input}.fai | "
        "awk '{{print $1\"\t0\t\"$2}}' >  {output}; "
        "rm {input}.fai"
