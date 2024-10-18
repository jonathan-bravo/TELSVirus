rule get_host_ref:
    output:
        f"{REF_DOWNLOADS}/{HOST_DOWNLOAD}",
    params:
        organism=HOST,
    conda:
        "../envs/get_host.yaml"
    benchmark:
        f"{BENCHDIR}/get_host_ref.benchmark"
    log:
        f"{LOGDIR}/get_host_ref_snakemake.log",
    shell:
        "python workflow/scripts/get_ref.py "
        "--organism '{params.organism}' "
        "--outfile {output}"
