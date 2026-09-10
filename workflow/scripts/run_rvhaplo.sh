#!/bin/bash
set -euo pipefail

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 OUTDIR REF_DIR SAM_DIR SAMPLE THREADS MAX_ALIGNMENTS MAX_ALIGNMENT_BASES" >&2
    exit 1
fi
outdir=$1
vir_indir=$2
sam_indir=$3
barcode=$4
threads=$5
max_alignments=$6
max_alignment_bases=$7

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Compatibility fixes are built into the pinned TELSVirus RVHaplo fork.
expected_revision=$(git -C "${REPO_ROOT}" rev-parse :RVHaplo)
actual_revision=$(git -C "${REPO_ROOT}/RVHaplo" rev-parse HEAD)
if [ "${actual_revision}" != "${expected_revision}" ]; then
    echo "ERROR: RVHaplo checkout does not match TELSVirus. Run git submodule sync --recursive and git submodule update --init --recursive." >&2
    exit 1
fi

mkdir -p "${outdir}"
abs_out="$(cd "${outdir}" && pwd)"
report="${abs_out}/rvhaplo_preflight.tsv"
printf 'sample\treference\talignments\treference_length\talignment_bases\tmax_alignments\tmax_alignment_bases\tdecision\treason\tsubgraphs\tthreads\n' > "${report}"

shopt -s nullglob
references=("${vir_indir%/}/"*.fasta)
if [ "${#references[@]}" -eq 0 ]; then
    echo "ERROR: No reference FASTA files found in ${vir_indir}" >&2
    exit 1
fi

for r in "${references[@]}"; do
    file=$(basename -- "$r")
    vir=${file%.fasta}
    sam_file="${sam_indir%/}/${barcode}_aligned_${vir}.sam"
    target_out="${abs_out}/rvhaplo_${barcode}_${vir}"
    mkdir -p "${target_out}"
    # Until this attempt succeeds, do not parse an old result from this directory.
    touch "${target_out}/rvhaplo-failed.flag"
    metrics=$(python "${SCRIPT_DIR}/rvhaplo_preflight.py" \
        --sam "${sam_file}" --reference "$r" \
        --max-alignments "${max_alignments}" \
        --max-alignment-bases "${max_alignment_bases}")
    IFS=$'\t' read -r read_count ref_length product decision reason <<< "${metrics}"

    sub_graph=1
    if [ "${read_count}" -gt 50000 ]; then
        sub_graph=$((read_count / 25000))
    fi
    if [ "${decision}" = skip ]; then
        sub_graph=0
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${barcode}" "${vir}" "${read_count}" "${ref_length}" "${product}" \
        "${max_alignments}" "${max_alignment_bases}" "${decision}" "${reason}" \
        "${sub_graph}" "${threads}" >> "${report}"
    echo "RVHaplo target=${vir} alignments=${read_count} reference_length=${ref_length} alignment_bases=${product} decision=${decision} reason=${reason} subgraphs=${sub_graph} threads=${threads}"

    if [ "${decision}" = skip ]; then
        printf '%s\n' "${reason}" > "${target_out}/rvhaplo-skipped.flag"
        rm -f "${target_out}/rvhaplo-failed.flag"
        continue
    fi
    rm -f "${target_out}/rvhaplo-skipped.flag" "${target_out}/viral-ref-too-long.flag"
    abs_sam=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${sam_file}")
    abs_ref=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$r")
    # Avoid accepting stale haplotypes if this attempt fails to produce output.
    rm -f "${target_out}/rvhaplo_haplotypes.fasta"
    if (cd "${REPO_ROOT}/RVHaplo" && bash rvhaplo.sh \
        -i "${abs_sam}" -r "${abs_ref}" -t "${threads}" -mq 0 \
        -sg "${sub_graph}" -o "${target_out}/"); then
        rm -f "${target_out}/rvhaplo-failed.flag"
    else
        echo "ERROR: RVHaplo failed for ${vir}; see console output." >&2
        exit 1
    fi
done
