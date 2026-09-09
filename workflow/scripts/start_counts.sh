input=$1
output=$2

echo $(gzip -dc "${input}" | wc -l)/4 | bc > "${output}"