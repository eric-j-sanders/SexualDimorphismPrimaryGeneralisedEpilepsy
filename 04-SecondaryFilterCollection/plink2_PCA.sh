#!/bin/bash

genotype=$1
output=$2

export PATH="$PATH:/hpf/largeprojects/struglis/eric/Programs/miniconda3/envs/gwas2024/bin:/hpf/largeprojects/struglis/eric/Programs/miniconda3/condabin:/usr/local/bin:/usr/local/sbin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/home/eric/bin:/hpf/largeprojects/struglis/eric/Programs/miniconda3/envs/gwas2024/bin:/hpf/largeprojects/struglis/eric/Programs/miniconda3/condabin"

num_samples=$(wc -l < "${genotype}.psam")
num_variants=$(wc -l < "${genotype}.pvar") 

if [[ $num_samples -lt $num_variants ]]; then
    n=$((num_samples - 1))
else
    n=$((num_variants - 1))
fi

plink2 --pfile $genotype --pca $n --out $output