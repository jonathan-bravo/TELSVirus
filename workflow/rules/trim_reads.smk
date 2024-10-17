rule trim_reads: # want to trim nanopore (8) + UMI (5) + illumina adaptors (24) 
    input:
       f"{OUTDIR}/{{sample}}.concat.fastq.gz"
    output:
        logfile = f"{OUTDIR}/{{sample}}.trimmed.log",
        trimmed_reads = f"{OUTDIR}/{{sample}}.trimmed.fastq.gz"
    params:
        crop = config["crop_len"],
        barcodes = config["barcodes"]
    conda:
        "../envs/trimming.yaml"
    log:
        f"{LOGDIR}/{{sample}}_trim_reads_snakemake.log"
    shell:
        "python workflow/scripts/trim.py "
        "--infile {input} "
        "--logfile {output.logfile} "
        "--outfile {output.trimmed_reads} "
        "--barcodes {params.barcodes} "
        "--crop {params.crop}"