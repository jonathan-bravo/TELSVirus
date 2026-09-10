rule gen_strain_db:
    input:
        VIRUSES,
    output:
        f"{REF_DOWNLOADS}/strain_db.tsv",
    log:
        f"{LOGDIR}/gen_strain_db_snakemake.log",
    benchmark:
        f"{BENCHDIR}/gen_strain_db.benchmark"
    conda:
        "../envs/strain_db.yaml"
    params:
        email=config["email"],
    shell:
        "python workflow/scripts/gen_strain_source_db.py "
        "--infile {input:q} "
        "--email {params.email:q} "
        "--outfile {output:q}"
