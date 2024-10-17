rule gen_strain_db:
    input:
        VIRUSES,
    output:
        "resources/downloaded_ref_files/strain_db.tsv",
    params:
        email=config["email"],
    conda:
        "../envs/strain_db.yaml"
    benchmark:
        f"{BENCHDIR}/gen_strain_db.benchmark"
    log:
        f"{LOGDIR}/gen_strain_db_snakemake.log",
    shell:
        "python workflow/scripts/gen_strain_source_db.py "
        "--infile {input} "
        "--email {params.email} "
        "--outfile {output}"
