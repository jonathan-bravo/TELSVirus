rule align_to_viruses_with_secondary:
    input:
        viruses=VIRUSES,
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
    output:
        temp(f"{OUTDIR}/{{sample}}_viruses.sam"),
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_ato_viruses.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_ato_viruses_snakemake.log",
    shell:
        "minimap2 "
        "--q-occ-frac 0 "
        "-f 0 "
        "-t {threads} "
        "-o {output} "
        "-a {input.viruses} {input.samples}"


rule viruses_sam_to_bam:
    input:
        f"{OUTDIR}/{{sample}}_viruses.sam",
    output:
        bam=f"{OUTDIR}/{{sample}}_viruses_sorted_sftclp.bam",
        idx=temp(f"{OUTDIR}/{{sample}}_reads_per_strain_samtools.idxstats"),
    params:
        cutoff=CUTOFF,
        tmp_bam=f"{OUTDIR}/{{sample}}_viruses_sorted.bam",
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_viral_sam2bam.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_viral_sam2bam_snakemake.log",
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {params.tmp_bam}; "
        "samtools index -@ {threads} {params.tmp_bam}; "
        "python workflow/scripts/sftclipper.py "
        "--cutoff {params.cutoff} "
        "--bam {params.tmp_bam} "
        "--outfile {output.bam}; "
        "samtools index -@ {threads} {output.bam}; "
        "samtools view -h -F 2048 -b {output.bam} | "
        "samtools idxstats - > {output.idx}; "
        "rm {params.tmp_bam} {params.tmp_bam}.bai"
