#!/bin/bash
# Make sure extended_ld_regions.txt is present locally

dataset=$1
prefix=$2

echo "0. Identify duplicate IDs and remove"

plink2 --pfile ${dataset} --list-duplicate-vars ids-only --out ${prefix}-dupcheck
if [ -s ${prefix}-dupcheck.dupvar ]; then
    cp ${prefix}-dupcheck.dupvar dups_to_rm.txt
else
    > dups_to_rm.txt
fi

echo "1. Listing non-autosomal chromosome SNPs, indels and ambiguous SNPs to delete"
awk '$1 !~ /^#/ && ($1 < 1 || $1 > 22 || $1 == "X")' ${dataset}.pvar |cut -f3 >snps_to_rm.txt
awk '$4 == "I" || $5 == "I" || $4 == "D" || $5 == "D"' ${dataset}.pvar |cut -f3 >>snps_to_rm.txt
awk '($4 == "A" && $5 == "T") || ($4 == "T" && $5 == "A") || ($4 == "G" && $5 == "C") || ($4 == "C" && $5 == "G")' ${dataset}.pvar |cut -f3 >>snps_to_rm.txt
awk '($4 != "A" && $4 != "T" && $4 != "G" && $4 != "C") || ($5 != "A" && $5 != "T" && $5 != "G" && $5 != "C")' ${dataset}.pvar |cut -f3 >>snps_to_rm.txt

cat dups_to_rm.txt >> snps_to_rm.txt

plink2 --pfile ${dataset} --exclude snps_to_rm.txt --make-pgen --out ${prefix}-bad_snps_rm --memory 10000


echo "2. Remove long LD region, rare alleles, and prune it"
if [ ! -e extended_ld_regions.txt ]; then ln -s ~/scripts/extended_ld_regions.txt; fi
plink2 --pfile ${prefix}-bad_snps_rm --exclude range extended_ld_regions.txt --make-pgen --out ${prefix}-longregionRemoved --memory 10000
#prune
plink2 --pfile ${prefix}-longregionRemoved --maf 0.05 --hwe 1e-4 0.001 keep-fewhet --set-all-var-ids @:# --rm-dup force-first --make-pgen --out ${prefix}-MAFlongregionRemoved --memory 10000
plink2 --pfile ${prefix}-MAFlongregionRemoved --indep-pairwise 1500 100 0.2 --out ${prefix}-pruned --memory 10000
plink2 --pfile ${prefix}-MAFlongregionRemoved --extract ${prefix}-pruned.prune.in --make-pgen --out ${prefix}-pruned --memory 10000

rm ${prefix}-bad_snps_rm* ${prefix}-longregionRemoved* ${prefix}-MAFlongregionRemoved*

