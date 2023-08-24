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
    sam_file=${sam_indir}${barcode}.aligned.${vir}.sam;
    read_count=$(samtools view -c -F 260 ${sam_file});
    if [ ${read_count} -gt 50000 ];
    then
        sub_graph=$(echo ${read_count}/25000 | bc)
    fi;
    scripts/RVHaplo/rvhaplo.sh \
    -i ${sam_file} \
    -r ${r} \
    -t ${threads} \
    -sg ${sub_graph} \
    -l 0 \
    -o ${outdir}rvhaplo_${barcode}_${vir}/;
done