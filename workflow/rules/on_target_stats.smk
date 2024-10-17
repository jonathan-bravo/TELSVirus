rule viruses_alignment_stats: # not counting supp alignments for `viral.stats`
    input:
        f"{OUTDIR}/{{sample}}.stats.viruses.sam"
    output:
        temp(f"{OUTDIR}/{{sample}}.all.viruses.samtools.idxstats")
    params:
        cutoff = config["sftclp_cutoff"],
        bam = f"{OUTDIR}/{{sample}}.stats.viruses.sorted.bam",
        sftbam = f"{OUTDIR}/{{sample}}.stats.viruses.sorted.sftclp.bam"
    conda:
        "../envs/alignment.yaml"
    threads:
        32
    benchmark:
        f"{BENCHDIR}/{{sample}}_stat_align.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_stat_align_snakemake.log"
    shell:
        "samtools view -h -F 2048 -Sb {input} | "
        "samtools sort -@ {threads} -o {params.bam}; "
        "samtools index -@ {threads} {params.bam}; "
        "python workflow/scripts/sftclipper.py "
        "--cutoff {params.cutoff} "
        "--bam {params.bam} "
        "--outfile {params.sftbam}; "
        "samtools index -@ {threads} {params.sftbam}; "
        "samtools idxstats {params.sftbam} > {output}; "
        "rm {params.bam} {params.bam}.bai; "
        "rm {params.sftbam} {params.sftbam}.bai"


rule merge_viral_alignment_stats:
    input:
        expand(f"{OUTDIR}/{{sample}}.all.viruses.samtools.idxstats", sample = SAMPLES)
    output:
        temp(f"{OUTDIR}/all.viruses.stats")
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/merge_viral_stats.benchmark"
    log:
        f"{LOGDIR}/merge_viral_stats_snakemake.log"
    shell:
        "python workflow/scripts/samtools_idxstats.py -i {input} -o {output}"


rule on_target_stats:
    input:
        viral_stats = f"{OUTDIR}/all.viruses.stats",
        host_stats = f"{OUTDIR}/host.removal.stats"
    output:
        f"{OUTDIR}/on.target.stats"
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/on_target_stats.benchmark"
    log:
        f"{LOGDIR}/on_target_stats_snakemake.log"
    shell:
        "python workflow/scripts/merge_stats.py "
        "--host_stats {input.host_stats} "
        "--viral_stats {input.viral_stats} "
        "--outfile {output}"