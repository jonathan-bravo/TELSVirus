rule pre_dedup_read_lengths:
    input:
        f"{OUTDIR}/{{sample}}.trimmed.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}.pre.dedup.rl.tsv",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_pre_dedup_rl.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_pre_dedup_rl_snakemake.log",
    shell:
        "python workflow/scripts/read_lengths.py "
        "--infile {input} "
        "--outfile {output}"


rule post_dedup_read_lengths:
    input:
        f"{OUTDIR}/{{sample}}.dedup.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}.post.dedup.rl.tsv",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_post_dedup_rl.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_post_dedup_rl_snakemake.log",
    shell:
        "python workflow/scripts/read_lengths.py "
        "--infile {input} "
        "--outfile {output}"
