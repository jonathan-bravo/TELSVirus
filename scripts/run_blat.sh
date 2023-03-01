outdir=$1
cluster_indir=$2

mkdir -p ${outdir};

for f in ${cluster_indir}*;
do
    file="$(basename -- $f)";
    cluster=${file%.rl.clusters.fasta.gz};
    blat ${f} ${f} ${outdir}${cluster}.psl
done