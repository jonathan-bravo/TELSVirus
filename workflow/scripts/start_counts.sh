input=$1
output=$2

echo $(gzcat ${input} | wc -l)/4 | bc > ${output}