#!/bin/bash
#SBATCH --job-name=merge_bplink     # Job name
#SBATCH --output=job_output/merge_bplink_%j.out  # Output log file (%j expands to job ID)
#SBATCH --error=job_output/merge_bplink_%j.err   # Error log file
#SBATCH --time=96:00:00      
#SBATCH --mem=300G               
#SBATCH --cpus-per-task=4       
# Given two pgen PLINK files with independent sample sets, merges them at the common SNPs only -- after checking for strand flip and/or ref/alt swapping
# This version deletes the intermediary files and outputs the merged "output_filename"
# I built this version specifically for two pgen datasets that already have matching psam column names -- #FID, IID, SEX, STATUS
# It is also expected that pvar files are formatted the same, and refers to the HG38 fasta file located in /hpf/largeprojects/struglis/eric
# Syntax: bash merge_bplink_files_v1_for_plink2.sh <dataset1_filename> <dataset2_filename> <output_filename> <temp_identifier>
# Example: bash merge_bplink_files_v1_for_plink2.sh 03_spit_for_sci_ctrls 38_ceu_cts_IIDupdated_18nov2014update 04_UK_SS_merged chr11
# All formats should be PLINK pgen format

#set -e

source /hpf/largeprojects/struglis/eric/Programs/miniconda3/bin/activate gwas2024

outname=$3
id=$4

echo "0. Replacing matching SNP IDs with rsid in Dataset 2"

awk '
	NR==FNR {
		if ($1 ~ /^#/) next         
		key = $1 ":" $2             
		idmap[key] = $3
		next
	}

	$1 ~ /^#/ { print; next }

	{
		key = $1 ":" $2
		if (key in idmap) $3 = idmap[key]
		print
	}
' ${1}.pvar ${2}.pvar > temp_${id}_step0_rsid2.pvar
cp ${2}.pgen temp_${id}_step0_rsid2.pgen
cp ${2}.psam temp_${id}_step0_rsid2.psam

echo "1. Removing Duplicate SNPs in Dataset1"

plink2 --pfile $1 --make-pgen --out temp_${id}_step1_nodup1 --snps-only --allow-extra-chr --rm-dup exclude-all

echo "2. Removing Duplicate SNPs in Dataset2"

plink2 --pfile temp_${id}_step0_rsid2 --make-pgen --out temp_${id}_step2_nodup2 --snps-only --allow-extra-chr --rm-dup exclude-all

echo "3. Listing Dataset1 SNPs"
cut -f3 temp_${id}_step1_nodup1.pvar > temp_${id}_step3.dataset_snplist.txt

echo "4. Extracting Dataset common SNPs from Dataset2"
plink2 --pfile temp_${id}_step2_nodup2 --extract temp_${id}_step3.dataset_snplist.txt --make-pgen --out temp_${id}_step4_dataset2_common_snps --memory 16000

echo "5. Listing Dataset2 common SNPs"
cut -f3 temp_${id}_step4_dataset2_common_snps.pvar > temp_${id}_step5.dataset2_common_snps.txt

echo "6. Extracting common SNPs from Dataset1"
plink2 --pfile temp_${id}_step1_nodup1 --extract temp_${id}_step5.dataset2_common_snps.txt --make-pgen --out temp_${id}_step6_dataset1_common_snps --memory 16000

echo "7. Convert to vcf"

# vcf will forcefully combine FID and IID via underscores
# remove future headache by replacing pre-existing underscores with hyphens
sed -i 's/_/-/g' temp_${id}_step6_dataset1_common_snps.psam
sed -i 's/_/-/g' temp_${id}_step4_dataset2_common_snps.psam

# --export is supposed to have "bgz" option, but it was producing corrupt .vcf.gz files
plink2 --pfile temp_${id}_step6_dataset1_common_snps --export vcf vcf-dosage=DS-force --out temp_${id}_step7_formerge1
plink2 --pfile temp_${id}_step4_dataset2_common_snps --export vcf vcf-dosage=DS-force --out temp_${id}_step7_formerge2

# Manually zip after producing
bgzip temp_${id}_step7_formerge1.vcf
bgzip temp_${id}_step7_formerge2.vcf

bcftools index temp_${id}_step7_formerge1.vcf.gz
bcftools index temp_${id}_step7_formerge2.vcf.gz

echo "8. Normalizing to reference fasta to hopefully allow for smooth merge"

{
  for i in {1..22}; do
    echo -e "${i}\tchr${i}";
  done
  echo -e "X\tchrX"
  echo -e "Y\tchrY"
  echo -e "MT\tchrM"
} > temp_${id}_step8_rename_chrs.txt

bcftools annotate --rename-chrs temp_${id}_step8_rename_chrs.txt -Oz -o temp_${id}_step8_formerge1_renamed.vcf.gz temp_${id}_step7_formerge1.vcf.gz
bcftools index -t temp_${id}_step8_formerge1_renamed.vcf.gz

bcftools annotate --rename-chrs temp_${id}_step8_rename_chrs.txt -Oz -o temp_${id}_step8_formerge2_renamed.vcf.gz temp_${id}_step7_formerge2.vcf.gz
bcftools index -t temp_${id}_step8_formerge2_renamed.vcf.gz

# Comparing to reference fasta seems broken

#bcftools +fixref temp_${id}_step8_formerge1_renamed.vcf.gz -- -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip -d | \
#bcftools norm -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m -both -Oz -o temp_${id}_step8_formerge1_normalized.vcf.gz
#bcftools sort temp_${id}_step8_formerge1_normalized.vcf.gz -Oz -o temp_${id}_step8_formerge1_normalized.sorted.vcf.gz
#bcftools index temp_${id}_step8_formerge1_normalized.sorted.vcf.gz

#bcftools +fixref temp_${id}_step8_formerge2_renamed.vcf.gz -- -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip -d | \
#bcftools norm -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m -both -Oz -o temp_${id}_step8_formerge2_normalized.vcf.gz
#bcftools sort temp_${id}_step8_formerge2_normalized.vcf.gz -Oz -o temp_${id}_step8_formerge2_normalized.sorted.vcf.gz
#bcftools index temp_${id}_step8_formerge2_normalized.sorted.vcf.gz

# Old broken way to reference fasta
#bcftools +fixref temp_${id}_step8_formerge1_renamed.vcf.gz -Oz -o temp_${id}_step8_formerge1_prepared.vcf.gz \
#	-- -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip -m swap -d
#bcftools +fixref temp_${id}_step8_formerge2_renamed.vcf.gz -Oz -o temp_${id}_step8_formerge2_prepared.vcf.gz \
#	-- -f /hpf/largeprojects/struglis/eric/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip -m swap -d
#bcftools index -t temp_${id}_step8_formerge1_prepared.vcf.gz
#bcftools index -t temp_${id}_step8_formerge2_prepared.vcf.gz

echo "9. Do merge using bcftools and convert back to plink"

bcftools isec -Oz -p temp_${id}_step8_int -c none -n=2 temp_${id}_step8_formerge1_renamed.vcf.gz temp_${id}_step8_formerge2_renamed.vcf.gz

bcftools merge -Oz -o temp_${id}_step9_merged.vcf.gz temp_${id}_step8_int/0000.vcf.gz temp_${id}_step8_int/0001.vcf.gz

echo "Preparing sex/status"

# Take original two psam files, get new file with just new-format ID's and a sex column
awk 'NR==1 {print "#IID SEX"; next} {O1=$1"_"$2; print O1, $3}' temp_${id}_step4_dataset2_common_snps.psam > temp_${id}_step9_finalsex.txt
awk 'NR>1 {O1=$1"_"$2; print O1, $3}' temp_${id}_step6_dataset1_common_snps.psam >> temp_${id}_step9_finalsex.txt

# Take original two psam files, get new file with just new-format ID's and a status column
awk 'NR==1 {print "#IID STATUS"; next} {O1=$1"_"$2; print O1, 2}' temp_${id}_step4_dataset2_common_snps.psam > temp_${id}_step9_finalstatus.txt
awk 'NR>1 {O1=$1"_"$2; print O1, 1}' temp_${id}_step6_dataset1_common_snps.psam >> temp_${id}_step9_finalstatus.txt

echo "Converting to plink"

# Convert merged vcf back to plink, using the prepared sex dataset and status dataset.
plink2 --vcf temp_${id}_step9_merged.vcf.gz dosage=DS --update-sex temp_${id}_step9_finalsex.txt --pheno temp_${id}_step9_finalstatus.txt --make-pgen --out temp_${id}_step9_merged

echo "10. Cleanup"

if [ -f "temp_${id}_step9_merged.pgen" ]; then
	mv temp_${id}_step9_merged.pgen ${outname}.pgen
	mv temp_${id}_step9_merged.pvar ${outname}.pvar
	mv temp_${id}_step9_merged.psam ${outname}.psam
	mv temp_${id}_step9_merged.log ${outname}.log
	rm -r temp_${id}*
	echo "Done."
else
	echo "There was a problem -- Leaving temp files for debugging."
fi
