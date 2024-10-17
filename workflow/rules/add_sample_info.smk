rule add_sample_info:
    input:
        f"{OUTDIR}/{{sample}}.pre.dedup.rl.tsv",
        f"{OUTDIR}/{{sample}}.post.dedup.rl.tsv",
        f"{OUTDIR}/{{sample}}.reads.per.strain.tsv",
        f"{OUTDIR}/{{sample}}.reads.per.strain.filtered.tsv",
        f"{OUTDIR}/{{sample}}.viral.targets.log",
        f"{OUTDIR}/{{sample}}.selected.viral.targets.log"
    output:
        touch(f"{OUTDIR}/{{sample}}_add_sample_info.done")
    params:
        run_id = config["run_id"],
        sample_id = "{sample}"
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_add_sample_info.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_sample_info_snakemake.log"
    shell:
        "python workflow/scripts/add_run_info.py "
        "--sample_id {params.sample_id} "
        "--run_id {params.run_id} "
        "--infiles {input}"