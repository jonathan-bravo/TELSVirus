mkdir $2

awk -v out=$2 'BEGIN {n=0;} /^>/ {
    if(n%1==0){file=sprintf(out"/chunk%d.fasta",n);}
    print >> file;
    n++;
    next;
} { print >> file; }' < $1


for f in $2/chunk*;
do
    line=$(head -n 1 $f | cut -d':' -f1);
    name=${line:1}.fasta;
    mv $f $2/$name;
done