git clone $1;
chmod +x RVHaplo/rvhaplo.sh;
chmod +x RVHaplo/src/*;
sed -i 's^./src/^scripts/RVHaplo/src/^' RVHaplo/rvhaplo.sh;
sed -i 's^./src/^scripts/RVHaplo/src/^' RVHaplo/src/out_haplotypes.py;
mv RVHaplo/ scripts/