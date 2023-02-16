configfile: "config.yaml"

READS = config["reads"]
OUTDIR = config["output"]
HOST_FILE = config["host_genome"]
BAITS = config["baits"]
VIRUSES = config["viral_genomes"]
BARCODES = [f for f in os.listdir(READS) if not f.startswith('.')]

all_input = [
    OUTDIR + "host.removal.stats",
    expand(OUTDIR + "{barcode}/{barcode}.mpileup", barcode = BARCODES)
]

rule all:
    input:
        all_input

rule concat_parts:
    input:
        READS + "{barcode}"
    output:
        temp(OUTDIR + "{barcode}/concat_{barcode}.fastq.gz")
    shell:
        "cat {input}/* > {output}"

## REMOVE HOST DNA
rule align_reads_to_host: # check minimap flags from telseq
    input:
        host = HOST_FILE,
        barcodes = OUTDIR + "{barcode}/concat_{barcode}.fastq.gz"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.host.sam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        12
    shell:
        "minimap2 -t {threads} -a {input.host} {input.barcodes} -o {output}"

rule host_sam_to_bam:
    input:
        OUTDIR + "{barcode}/{barcode}.host.sam"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.host.sorted.bam")
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
        OUTDIR + "{barcode}/{barcode}.host.sorted.bam"
    output:
        idx = temp(OUTDIR + "{barcode}/{barcode}.remove.host.samtools.idxstats"),
        bam = temp(OUTDIR + "{barcode}/{barcode}.host.removed.sorted.bam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools index {input} && "
        "samtools idxstats {input} > {output.idx}; "
        "samtools view -h -f 4 -b {input} -o {output.bam}; "
        "rm {input}.bai"

rule host_removal_stats: ## ADD THE SCRIPT FROM AMR++, maybe modify
    input:
        expand(OUTDIR + "{barcode}/{barcode}.remove.host.samtools.idxstats", barcode = BARCODES)
    output:
        OUTDIR + "host.removal.stats"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/samtools_idxstats.py -i {input} -o {output}"

rule non_host_reads:
    input:
        OUTDIR + "{barcode}/{barcode}.host.removed.sorted.bam"
    output:
        OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bedtools/2.30.0"
    shell:
        "bedtools bamtofastq -i {input} -fq {output}"

## CALCULATE GENOME FRACTION
rule align_to_baits:
    input:
        baits = VIRUSES,
        barcodes = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.baits.sam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        12
    shell:
        "minimap2 -t {threads} -a {input.baits} {input.barcodes} -o {output}"

rule baits_sam_to_bam:
    input:
        OUTDIR + "{barcode}/{barcode}.baits.sam"
    output:
        OUTDIR + "{barcode}/{barcode}.baits.sorted.bam"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}; "
        "samtools index {output}"

rule temp_pileup:
    input:
        bam = OUTDIR + "{barcode}/{barcode}.baits.sorted.bam",
        baits = BAITS
    output:
        OUTDIR + "{barcode}/{barcode}.mpileup"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools mpileup -aa --output-QNAME -o {output} {input.bam}"



# -f {input.baits}