configfile: "config.yaml"

READS = config["reads"]
OUTDIR = config["output"]
HOST_FILE = config["host_genome"]
VIRUSES = config["viral_genomes"]
BARCODES = [f for f in os.listdir(READS) if not f.startswith('.')]

all_input = [
    OUTDIR + "host.removal.stats",
    expand(OUTDIR + "{barcode}/{barcode}.reads.per.strain.samtools.idxstats", barcode = BARCODES),
    expand(OUTDIR + "{barcode}/rvhaplo.done", barcode = BARCODES)
    ##expand(OUTDIR + "{barcode}/{barcode}_strainline_out/", barcode = BARCODES)
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

# rule get_strainline:
#     output:
#         touch("get_strainline.done")
#     params:
#         repo = config["strainline_repo"]
#     conda:
#         "envs/git.yaml"
#     envmodules:
#         "git/2.30.1"
#     shell:
#         "scripts/get_strainline.sh {params.repo}; "
#         "scripts/get_daccord.sh" 

# rule link_daccord:
#     input:
#         "get_strainline.done"
#     output:
#         touch("daccord_linked.done")
#     conda:
#         "envs/strainline.yaml"
#     shell:
#         "ln -fs scripts/daccord/bin/daccord "
#         "$CONDA_PREFIX/bin/daccord"

rule concat_parts:
    input:
        READS + "{barcode}"
    output:
        temp(OUTDIR + "{barcode}/concat_{barcode}.fastq.gz")
    shell:
        "cat {input}/* > {output}"

## REMOVE ADAPTERS
## FIRST 5bp NEED TO BE TRIMMED OFF
## TRIMMOMATIC

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
        "minimap2 --secondary=no "
        "-t {threads} "
        "-o {output} "
        "-a {input.host} {input.barcodes}"

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
        temp(OUTDIR + "host.removal.stats")
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
        temp(OUTDIR + "{barcode}/{barcode}.viruses.sam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        32
    shell:
        "minimap2 "
        "-t {threads} "
        "-o {output} "
        "-a {input.viruses} {input.barcodes}"

rule align_to_viruses_for_stats:
    input:
        viruses = VIRUSES,
        barcodes = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.stats.viruses.sam")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        32
    shell:
        "minimap2 "
        "--secondary=no "
        "-t {threads} "
        "-o {output} "
        "-a {input.viruses} {input.barcodes}"

rule viruses_sam_to_bam:
    input:
        OUTDIR + "{barcode}/{barcode}.viruses.sam"
    output:
        bam = OUTDIR + "{barcode}/{barcode}.viruses.sorted.bam",
        idx = temp(OUTDIR + "{barcode}/{barcode}.reads.per.strain.samtools.idxstats")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        32
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output.bam}; "
        "samtools index {output.bam}; "
        "samtools idxstats {output.bam} > {output.idx}"

rule reads_per_strain:
    input:
        OUTDIR + "{barcode}/{barcode}.reads.per.strain.samtools.idxstats"
    output:
        OUTDIR + "{barcode}/{barcode}.reads.per.strain.tsv"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/summarize_idxstats.py "
        "--infile {input} "
        "--outfile {output}"

rule viruses_alignment_stats:
    input:
        OUTDIR + "{barcode}/{barcode}.stats.viruses.sam"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.all.viruses.samtools.idxstats")
    params:
        bam = "{barcode}.stats.viruses.sorted.bam"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    threads:
        32
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {params.bam}; "
        "samtools index {params.bam}; "
        "samtools idxstats {params.bam} > {output}; "
        "rm {params.bam}.bai; "
        "rm {params.bam}"

rule merge_alignment_stats: ## ADD THE SCRIPT FROM AMR++, maybe modify
    input:
        expand(OUTDIR + "{barcode}/{barcode}.all.viruses.samtools.idxstats", barcode = BARCODES)
    output:
        temp(OUTDIR + "all.viruses.stats")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/samtools_idxstats.py -i {input} -o {output}"

rule on_target_stats:
    input:
        viral_stats = OUTDIR + "all.viruses.stats",
        host_stats = OUTDIR + "host.removal.stats"
    output:
        OUTDIR + "on.target.stats"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/merge_stats.py "
        "--host_stats {input.host_stats}"
        "--viral_stats {input.viral_stats}"
        "--outfile {output}"

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
        "samtools mpileup --output-QNAME {input} | "
        "awk '{{print $1\"\t\"$4}}' - > {output}"

rule all_virus_bed:
    input:
        VIRUSES
    output:
        temp(OUTDIR + "all.viral.targets.bed")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.9"
    shell:
        "samtools faidx {input}; "
        "cat {input}.fai | "
        "awk '{{print $1\"\t0\t\"$2}}' >  {output}; "
        "rm {input}.fai"

rule gen_strain_db:
    input:
        VIRUSES
    output:
        OUTDIR + "strain_db.json"
    params:
        email = config["email"]
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/gen_strain_db.py "
        "--infile {input} "
        "--email {params.email} "
        "--outfile {output}"

rule find_viral_tagets:
    input:
        pileup = OUTDIR + "{barcode}/{barcode}.mpileup",
        all_viruses_bed = OUTDIR + "all.viral.targets.bed",
        strain_db = OUTDIR + "strain_db.json"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.viral.targets.bed")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8"
    shell:
        "scripts/find_viral_targets.py "
        "--mpileup {input.pileup} "
        "--bed {input.all_viruses_bed} "
        "--outfile {output}"
        
        # "grep -Ff {params.target_list} - > {output}; "

rule get_viral_genomes:
    input:
        fasta = VIRUSES,
        bed = OUTDIR + "{barcode}/{barcode}.viral.targets.bed"
    output:
        temp(OUTDIR + "{barcode}/{barcode}.viral.target.genomes.fasta")
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bedtools/2.30.0"
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output}"

rule split_viral_genomes:
    input:
        OUTDIR + "{barcode}/{barcode}.viral.target.genomes.fasta"
    output:
        touch(OUTDIR + "{barcode}/viral_refs.done")
    params:
        outdir = OUTDIR + "{barcode}/viral_refs/"
    shell:
        "scripts/split_target_viruses.sh {input} {params.outådir}"

# [x] align non.host.fastas to ALL the split viral genomes
# [ ] Then run rvhaplo on THOSE sam files
# [x] can then mark the other sam file as temp()
# [ ] minimap2 -a -x map-ont MT269879.fasta cat_barcode04.fastq.gz -o barcode04_2_map2.sam

rule align_to_target_virus:
    input:
        OUTDIR + "{barcode}/viral_refs.done",
        reads = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
    output:
        touch(OUTDIR + "{barcode}/align_to_targets.done")
    params:
        indir = OUTDIR + "{barcode}/viral_refs/",
        outdir = OUTDIR + "{barcode}/target_aligned/",
        barcode = "{barcode}"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "minimap2/2.24"
    threads:
        32
    shell:
        "scripts/align_to_targets.sh "
        "{params.outdir} "
        "{params.indir} "
        "{params.barcode} "
        "{input.reads} "
        "{threads}"

rule run_rvhaplo:
    input:
        "get_rvhaplo.done",
        OUTDIR + "{barcode}/align_to_targets.done"
    output:
        touch(OUTDIR + "{barcode}/rvhaplo.done")
    params:
        barcode = "{barcode}",
        vir_indir = OUTDIR + "{barcode}/viral_refs/",
        sam_indir = OUTDIR + "{barcode}/target_aligned/",
        outdir = OUTDIR + "{barcode}/rvhaplo_out/"
    conda:
        "envs/rvhaplo.yaml"
    threads:
        10
    shell:
        "scripts/run_rvhaplo.sh "
        "{params.outdir} "
        "{params.vir_indir} "
        "{params.sam_indir} "
        "{params.barcode} "
        "{threads}"

# LAS error, not able to read fatsa file, need to figure out why
# On Hipergator the `daccord` command isn't being found in some samples
## even after the command has been linked correctly to the bin

# rule run_strainline: ## need to figure out how to link daccord to strainline conda bin...
#     input:
#         "daccord_linked.done",
#         fasta = OUTDIR + "{barcode}/{barcode}.non.host.fastq.gz"
#     output:
#         OUTDIR + "{barcode}/{barcode}_strainline_out/"
#     threads:
#         10
#     conda:
#         "envs/strainline.yaml"
#     shell:
#         "scripts/Strainline/src/strainline.sh "
#         "-i {input.fasta} "
#         "-o {output} "
#         "-t {threads} "
#         "-p ont"