rule reads_per_strain:
    input:
        stats=f"{OUTDIR}/{{sample}}_reads_per_strain_samtools.idxstats",
        strain_db=STRAIN_DB,
    output:
        f"{OUTDIR}/{{sample}}_reads_per_strain.tsv",
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_reads_per_strain.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_reads_per_strain_snakemake.log",
    shell:
        "python workflow/scripts/summarize_idxstats.py "
        "--infile {input.stats} "
        "--strains {input.strain_db} "
        "--outfile {output}"


rule filter_reads_per_strain:
    input:
        f"{OUTDIR}/{{sample}}_reads_per_strain.tsv",
    output:
        f"{OUTDIR}/{{sample}}_reads_per_strain_filtered.tsv",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_filtered_reads_per_strain.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_filtered_reads_per_strain_snakemake.log",
    shell:
        "cat {input} | "
        "awk -F'\t' '{{ if($3 > 0) print}}' > {output}"
