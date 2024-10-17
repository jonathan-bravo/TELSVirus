rule align_reads_to_host:
    input:
        host = HOST_FILE,
        samples = f"{OUTDIR}/{{sample}}.dedup.fastq.gz",
    output:
        temp(f"{OUTDIR}/{{sample}}.host.sam")
    conda:
        "../envs/alignment.yaml"
    threads:
        32
    log:
        f"{LOGDIR}/{{sample}}_ato_host_snakemake.log"
    shell:
        "minimap2 --secondary=no "
        "-t {threads} "
        "-o {output} "
        "-a {input.host} {input.samples}"


rule host_sam_to_bam:
    input:
        f"{OUTDIR}/{{sample}}.host.sam"
    output:
        temp(f"{OUTDIR}/{{sample}}.host.sorted.bam")
    conda:
        "../envs/alignment.yaml"
    threads:
        32
    log:
        f"{LOGDIR}/{{sample}}_host_sam2bam_snakemake.log"
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}"


rule remove_host_dna: # not counting supp alignments for `host.removal`
    input:
        f"{OUTDIR}/{{sample}}.host.sorted.bam"
    output:
        idx = temp(f"{OUTDIR}/{{sample}}.remove.host.samtools.idxstats"),
        bam = temp(f"{OUTDIR}/{{sample}}.host.removed.sorted.bam")
    conda:
        "../envs/alignment.yaml"
    threads:
        10
    log:
        f"{LOGDIR}/{{sample}}_host_removal_snakemake.log"
    shell:
        "samtools index {input} && "
        "samtools view -h -F 2048 -b {input} | "
        "samtools idxstats - > {output.idx}; "
        "samtools view -h -f 4 -b {input} -o {output.bam}; "
        "rm {input}.bai"


rule non_host_reads:
    input:
        f"{OUTDIR}/{{sample}}.host.removed.sorted.bam"
    output:
        f"{OUTDIR}/{{sample}}.non.host.fastq.gz"
    # params:
    #     fq = f"{OUTDIR}/{{sample}}.non.host.fastq"
    conda:
        "../envs/alignment.yaml"
    threads:
        4
    log:
        f"{LOGDIR}/{{sample}}_non_host_reads_snakemake.log"
    shell:
        # "samtools fastq {input} > {params.fq}; "
        # "bgzip -@ {threads} {params.fq}"
        "samtools fastq {input} | bgzip -@ {threads} -c > {output}"
