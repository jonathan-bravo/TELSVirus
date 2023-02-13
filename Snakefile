configfile: "config.yaml"

READS = config["reads"]
OUTDIR = config["output"]
HOST_FILE = config["host_genome"]
TAGS, = glob_wildcards(READS + "*_pass_barcode{tag}_*.fastq.gz")

all_input = []

rule all:
    input:
        all_input


rule concat_parts:
    input:
        READS + "*_pass_barcode{tag}_*.fastq.gz"
    output:
        OUTDIR + "barcode{tag}/concat_barcode{tag}.fastq.gz"
    shell:
        "cat {input} > concat_barcode{}"


# from AMR ++ will modify
## REMOVE HOST DNA
rule align_reads_to_host: # check minimap flags from telseq
    input:
        host = HOST_FILE,
        barcodes = OUTDIR + "barcode{tag}/concat_barcode{tag}.fastq.gz"
    output:
        temp(OUTDIR + "barcode{tag}/barcode{tag}.host.sam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        12
    shell:
        "minimap2 -t {threads} {input.host} {input.barcodes} -o {output}"


rule host_sam_to_bam:
    input:
        OUTDIR + "barcode{tag}/barcode{tag}.host.sam"
    output:
        temp(OUTDIR + "barcode{tag}/barcode{tag}.host.sorted.bam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}"


rule remove_host_dna:
    input:
        OUTDIR + "barcode{tag}/barcode{tag}.host.sorted.bam"
    output:
        idx = OUTDIR + "barcode{tag}/RemoveHostDNA/barcode{tag}.samtools.idxstats",
        bam = OUTDIR + "barcode{tag}/RemoveHostDNA/barcode{tag}.host.removed.sorted.bam"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools index {input} && "
        "samtools idxstats {input} > {output.idx}; "
        "samtools view -h -f 4 -b {input} -o {output.bam}"


rule host_removal_stats: ## ADD THE SCRIPT FROM AMR++, maybe modify
    input:
        expand(OUTDIR + "barcode{tag}/RemoveHostDNA/barcode{tag}.samtools.idxstats", tag = TAGS)
    output:
        OUTDIR + "barcode{tag}/RemoveHostDNA/host.removal.stats"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/samtools_idxstats.py -i {input} -o {output}"


rule non_host_reads:
    input:
        OUTDIR + "barcode{tag}/RemoveHostDNA/barcode{tag}.host.removed.sorted.bam"
    output:
        OUTDIR + "barcode{tag}/NonHostReads/barcode{tag}.non.host.fastq.gz",
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bedtools/2.30.0"
    shell:
        "bedtools bamtofastq -i {input} -fq {output}"
