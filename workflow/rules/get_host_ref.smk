rule get_host_ref:
    output:
        f"{REF_DOWNLOADS}/{HOST_DOWNLOAD}",
    log:
        f"{LOGDIR}/get_host_ref_snakemake.log",
    benchmark:
        f"{BENCHDIR}/get_host_ref.benchmark"
    conda:
        "../envs/get_host.yaml"
    params:
        organism=HOST,
    shell:
        "python workflow/scripts/get_ref.py "
        "--organism {params.organism:q} "
        "--outfile {output:q}"
