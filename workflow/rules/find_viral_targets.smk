rule mpileup:
    input:
        f"{OUTDIR}/{{sample}}.viruses.sorted.sftclp.bam"
    output:
        f"{OUTDIR}/{{sample}}.mpileup"
    conda:
        "../envs/alignment.yaml"
    threads:
        10
    log:
        f"{LOGDIR}/{{sample}}_mpileup_snakemake.log"
    shell:
        "samtools mpileup --output-QNAME {input} | "
        "awk '{{print $1\"\t\"$4}}' - > {output}"


rule find_viral_targets:
    input:
        pileup = f"{OUTDIR}/{{sample}}.mpileup",
        all_viruses_bed = f"{OUTDIR}/all.viral.targets.bed",
        strain_db = STRAIN_DB
    output:
        bed = temp(f"{OUTDIR}/{{sample}}.viral.targets.bed"),
        logfile = f"{OUTDIR}/{{sample}}.viral.targets.log",
        selected_log = f"{OUTDIR}/{{sample}}.selected.viral.targets.log"
    params:
        email = config["email"]
    conda:
        "../envs/alignment.yaml"
    log:
        f"{LOGDIR}/{{sample}}_viral_targets_snakemake.log"
    shell:
        "python workflow/scripts/find_viral_targets.py "
        "--mpileup {input.pileup} "
        "--bed {input.all_viruses_bed} "
        "--strains {input.strain_db} "
        "--logfile {output.logfile} "
        "--outfile {output.bed}; "
        "python workflow/scripts/add_viral_segments.py "
        "--infile {output.logfile} "
        "--email {params.email}; "
        "awk -F'\t' '{{ if($3 > 80.0) print}}' {output.logfile} > {output.selected_log}"


rule get_viral_genomes:
    input:
        fasta = VIRUSES,
        bed = f"{OUTDIR}/{{sample}}.viral.targets.bed"
    output:
        f"{OUTDIR}/{{sample}}.viral.target.genomes.fasta"
    params:
        succeed = f"{OUTDIR}/{{sample}}.VIRAL.TARGETS.FOUND",
        fail = f"{OUTDIR}/{{sample}}.NO.VIRAL.TARGETS"
    conda:
        "../envs/alignment.yaml"
    log:
        f"{LOGDIR}/{{sample}}_get_viral_genomes_snakemake.log"
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output}; "
        "if grep -q . {output}; then touch {params.succeed}; else touch {params.fail}; fi"