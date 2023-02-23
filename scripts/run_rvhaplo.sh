outdir=$1
vir_indir=$2
sam_indir=$3
barcode=$4
threads=$5

mkdir -p ${outdir};

for r in ${vir_indir}*;
do
    file="$(basename -- $r)";
    vir=${file%.fasta};
    scripts/RVHaplo/rvhaplo.sh \
    -i ${sam_indir}${barcode}.aligned.${vir}.sam \
    -r ${r} \
    -t ${threads} \
    -l 0 \
    -o ${outdir}rvhaplo_${barcode}_${vir}/;
done