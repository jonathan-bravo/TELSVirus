#!/bin/bash

outdir=$1
vir_indir=$2
sam_indir=$3
barcode=$4
threads=$5
sub_graph=1

mkdir -p ${outdir};

for r in ${vir_indir}*;
do
    file="$(basename -- $r)";
    vir=${file%.fasta};
    sam_file=${sam_indir}${barcode}_aligned_${vir}.sam;

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

    # Run RVHaplo, but don't fail if it errors
    bash RVHaplo/rvhaplo.sh \
        -i ${sam_file} \
        -r ${r} \
        -t ${threads} \
        -sg ${sub_graph} \
        -l 0 \
        -o ${outdir}rvhaplo_${barcode}_${vir}/ || echo "RVHaplo failed for ${vir}, continuing..."
done

# Always exit successfully - parse_rvhaplo.py will handle empty results
exit 0