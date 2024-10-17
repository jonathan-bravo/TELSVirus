rule reads_per_strain:
    input:
        stats = f"{OUTDIR}/{{sample}}.reads.per.strain.samtools.idxstats",
        strain_db = STRAIN_DB
    output:
        f"{OUTDIR}/{{sample}}.reads.per.strain.tsv"
    conda:
        "../envs/alignment.yaml"
    log:
        f"{LOGDIR}/{{sample}}_reads_per_strain_snakemake.log"
    shell:
        "python workflow/scripts/summarize_idxstats.py "
        "--infile {input.stats} "
        "--strains {input.strain_db} "
        "--outfile {output}"


rule filter_reads_per_strain:
    input:
        f"{OUTDIR}/{{sample}}.reads.per.strain.tsv"
    output:
        f"{OUTDIR}/{{sample}}.reads.per.strain.filtered.tsv"
    conda:
        "../envs/default.yaml"
    log:
        f"{LOGDIR}/{{sample}}_filtered_reads_per_strain_snakemake.log"
    shell:
        "cat {input} | "
        "awk -F'\t' '{{ if($3 > 0) print}}' > {output}"