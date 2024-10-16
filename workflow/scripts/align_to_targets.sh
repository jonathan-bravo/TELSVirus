outdir=$1
indir=$2
barcode=$3
reads=$4
threads=$5

mkdir -p ${outdir};

for r in ${indir}*;
do
    file="$(basename -- $r)";
    vir=${file%.fasta};
    minimap2 \
    --secondary=no \
    -t ${threads} \
    -a ${r} ${reads} \
    -o ${outdir}${barcode}_aligned_${vir}.sam;
done