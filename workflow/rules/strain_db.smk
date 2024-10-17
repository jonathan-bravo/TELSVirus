rule gen_strain_db:
    input:
        VIRUSES
    output:
        f"{OUTDIR}/downloaded_ref_files/strain_db.tsv"
    params:
        email = config["email"]
    conda:
        "../envs/strain_db.yaml"
    log:
        f"{LOGDIR}/gen_strain_db_snakemake.log"
    shell:
        "python workflow/scripts/gen_strain_source_db.py "
        "--infile {input} "
        "--email {params.email} "
        "--outfile {output}"
    