input=$1
output=$2

echo $(zcat ${input} | wc -l)/4 | bc > ${output}