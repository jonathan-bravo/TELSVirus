# #!/bin/bash

# outdir=$1
# vir_indir=$2
# sam_indir=$3
# barcode=$4
# threads=$5
# sub_graph=1

# mkdir -p ${outdir};

# for r in ${vir_indir}*;
# do
#     file="$(basename -- $r)";
#     vir=${file%.fasta};
#     sam_file=${sam_indir}${barcode}_aligned_${vir}.sam;

#     # Skip if SAM file doesn't exist or is empty
#     if [ ! -s "${sam_file}" ]; then
#         echo "Skipping ${vir}: SAM file not found or empty"
#         continue
#     fi

#     read_count=$(samtools view -c -F 260 ${sam_file} 2>/dev/null || echo "0");

#     # Skip if no aligned reads
#     if [ "${read_count}" -eq 0 ]; then
#         echo "Skipping ${vir}: no aligned reads"
#         continue
#     fi

#     if [ ${read_count} -gt 50000 ];
#     then
#         sub_graph=$(echo ${read_count}/25000 | bc)
#     fi;

#     # Run RVHaplo, but don't fail if it errors
#     bash RVHaplo/rvhaplo.sh \
#         -i ${sam_file} \
#         -r ${r} \
#         -t ${threads} \
#         -sg ${sub_graph} \
#         -l 0 \
#         -o ${outdir}rvhaplo_${barcode}_${vir}/ || echo "RVHaplo failed for ${vir}, continuing..."
# done

# # Always exit successfully - parse_rvhaplo.py will handle empty results
# exit 0

outdir=$1
vir_indir=$2
sam_indir=$3
barcode=$4
threads=$5
ref_length_limit=$6
sub_graph=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PATCH_FILE="${REPO_ROOT}/resources/rvhaplo-threadpool.patch"

# Apply patch to RVHaplo submodule if not already applied
(cd "${REPO_ROOT}/RVHaplo" && git apply --ignore-whitespace "${PATCH_FILE}" 2>/dev/null || true)

mkdir -p ${outdir};

for r in ${vir_indir}*;
do
    file="$(basename -- $r)";
    vir=${file%.fasta};
    sam_file=${sam_indir}${barcode}_aligned_${vir}.sam;
    ref_length=$(awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq $0} END {if (seq) print length(seq)}' $r)


    # Skip if SAM file doesn't exist or is empty
    if [ ! -s "${sam_file}" ]; then
        echo "Skipping ${vir}: SAM file not found or empty"
        continue
    fi

    read_count=$(samtools view -c -F 260 ${sam_file} 2>/dev/null || echo "0");

    # Skip if no aligned reads
    if [ "${read_count}" -eq 0 ]; then
        echo "Skipping ${vir}: no aligned reads"
        continue
    fi

    if [ ${read_count} -gt 50000 ];
    then
        sub_graph=$(echo ${read_count}/25000 | bc)
    fi;

    # Resolve to absolute paths before cd-ing into RVHaplo source directory
    abs_sam=$(readlink -f "${sam_file}")
    abs_ref=$(readlink -f "${r}")
    abs_out=$(readlink -f "${outdir}")

    if [ "$ref_length" -ge "$ref_length_limit" ]; then
        mkdir -p ${abs_out}/rvhaplo_${barcode}_${vir}/
        touch ${abs_out}/rvhaplo_${barcode}_${vir}/viral-ref-too-long.flag
        continue
    fi

    # Run RVHaplo from its source directory so ./src/ relative paths resolve correctly
    # old `-l 0` meaning that all values are included in clustering, default `-l 50`
    ## trying default, value but will have to tune parameter depending on data 
    (cd RVHaplo && bash rvhaplo.sh \
        -i "${abs_sam}" \
        -r "${abs_ref}" \
        -t ${threads} \
        -sg ${sub_graph} \
        -o "${abs_out}/rvhaplo_${barcode}_${vir}/") || echo "RVHaplo failed for ${vir}, continuing..."
done

# Always exit successfully - parse_rvhaplo.py will handle empty results
exit 0