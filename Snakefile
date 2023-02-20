configfile: "config.yaml"

READS = config["reads"]
OUTDIR = config["output"]
HOST_FILE = config["host_genome"]
VIRUSES = config["viral_genomes"]
BARCODES = [f for f in os.listdir(READS) if not f.startswith('.')]

all_input = [
    #"get_rvhaplo.done",
    #"get_strainline.done",
    #"daccord_linked.done",
    OUTDIR + "host.removal.stats",
    #expand(OUTDIR + "{barcode}/{barcode}.mpileup", barcode = BARCODES),
    expand(OUTDIR + "{barcode}/{barcode}_rvhaplo_out/", barcode = BARCODES),
    expand(OUTDIR + "{barcode}/{barcode}_strainline_out/", barcode = BARCODES)
]

rule all:
    input:
        all_input

rule get_rvhaplo:
    output:
        touch("get_rvhaplo.done")
    params:
        repo = config["rvhaplo_repo"]
    conda:
        "envs/git.yaml"
    envmodules:
        "git/2.30.1"
    shell:
        "scripts/get_rvhaplo.sh {params.repo}"

rule get_strainline:
    output:
        touch("get_strainline.done")
    params:
        repo = config["strainline_repo"]
    conda:
        "envs/git.yaml"
    envmodules:
        "git/2.30.1"
    shell:
        "scripts/get_strainline.sh {params.repo}; "
        "scripts/get_daccord.sh" 

rule link_daccord:
    input:
        "get_strainline.done"
    output:
        touch("daccord_linked.done")
    conda:
        "envs/strainline.yaml"
    shell:
        "ln -fs scripts/daccord/bin/daccord "
        "$CONDA_PREFIX/bin/daccord"

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
        32
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
        32
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
rule align_to_viruses:
    input:
        viruses = VIRUSES,
        barcodes = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    output:
        OUTDIR + "{barcode}/{barcode}.viruses.sam"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        32
    shell:
        "minimap2 -t {threads} -a {input.viruses} {input.barcodes} -o {output}"

rule viruses_sam_to_bam:
    input:
        OUTDIR + "{barcode}/{barcode}.viruses.sam"
    output:
        OUTDIR + "{barcode}/{barcode}.viruses.sorted.bam"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        32
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}; "
        "samtools index {output}"

rule mpileup:
    input:
        OUTDIR + "{barcode}/{barcode}.viruses.sorted.bam"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.mpileup")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        10
    shell:
        "samtools mpileup -d 1 --output-QNAME {input} | "
        "awk '{{print $1\"\t\"$4}}' - > {output}"

rule all_virus_bed:
    input:
        VIRUSES
    output:
        OUTDIR + "all.viral.targets.bed"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    shell:
        "samtools faidx {input}; "
        "cat {input}.fai | "
        "awk '{{print $1\"\t0\t\"$2}}' >  {output}; "
        "rm {input}.fai"

rule find_viral_tagets:
    input:
        pileup = OUTDIR + "{barcode}/{barcode}.mpileup",
        all_viruses_bed = OUTDIR + "all.viral.targets.bed"
    output:
        OUTDIR + "{barcode}/{barcode}.viral.targets.bed"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/find_viral_targets.py "
        "--mpileup {input.pileup} "
        "--bed {input.all_viruses_bed}"
        "--outfile {output}"
        
        # "grep -Ff {params.target_list} - > {output}; "

rule get_viral_genomes:
    input:
        fasta = VIRUSES,
        bed = OUTDIR + "{barcode}/{barcode}.viral.targets.bed"
    output:
        OUTDIR + "{barcode}/{barcode}.viral.target.genomes.fasta"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bedtools/2.30.0"
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output}"

rule run_rvhaplo:
    input:
        "get_rvhaplo.done",
        sam = OUTDIR + "{barcode}/{barcode}.host.sam",
        viral_ref = OUTDIR + "{barcode}/{barcode}.viral.target.genomes.fasta"
    output:
        OUTDIR + "{barcode}/{barcode}_rvhaplo_out/"
    conda:
        "envs/rvhaplo.yaml"
    shell:
        "scripts/RVHaplo/rvhaplo.sh "
        "-i {input.sam} "
        "-r {inpput.viral_ref}"
        "-l 0 "
        "-o {output}"

rule run_strainline: ## need to figure out how to link daccord to strainline conda bin...
    input:
        "daccord_linked.done",
        fasta = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    output:
        OUTDIR + "{barcode}/{barcode}_strainline_out/"
    threads:
        10
    conda:
        "envs/strainline.yaml"
    shell:
        "scripts/Strainline/src/strainline.sh "
        "-i {input.fasta} "
        "-o {output} "
        "-t {threads} "
        "-p ont"