# Purpose

When the codes in this folder are executed, the BIOJUME case data are filtered to European-ancestry individuals with the probable or definite JME phenotype assignment, and have some formatting changes to facilitate merging with controls.

# Input files

* A set of bed/bim/fam files formatted for plink are assumed to be located with the relative path prefix "rawdata/BIOJUME_R1234567", representing the genotyped (not imputed) data from the BIOJUME case population (before imputation). These data have already been put through the QC process as described in the paper.

* A set of 23 .vcf.gz files representing the imputed genotypes from the BIOJUME case population, is assumed to be located with the relative pathing "rawdata/BIOJUME_R1234567.rsq08.chr${chr}.vcf.gz"

* A phenotype file with columns corresponding to sample ID and covariates measured via patient questionaire or phenotyping committee (including phenotyping committee diagnostic assignments) is assumed to be located with the relative path "rawdata/BIOJUMEPal_DATA_2025-10-08_1639.csv".

# Primary outputs

plink2 format files for both the genotyped and imputed data, that have been subset to individuals with probable/definite JME:

* pgen/psam/pvar files with relative pathing prefix intdata/09-BIOJUME_prepared_genotype
* pgen/psam/pvar files with relative pathing prefix intdata/09-BIOJUME_prepared_imputed

# Code

## Subset and Format Genotype Data

Using R, create file listing probable/definite JME IDs:

```{r}
data = read.csv('rawdata/BIOJUMEPal_DATA_2025-10-08_1639.csv',na.strings=c("","N/A","na","NA"))
jme = data[data$phenotyping_classification %in% c(3,4),] # 3=probable, 4=definite JME

fam1 = read.table('rawdata/BIOJUME_R1234567.fam')

pheno = data.frame(FID=NA,study_id=jme[,'study_id'])
pheno$FID = fam1$V1[match(pheno$study_id, fam1$V2)]

pheno$female = fam1$V5[match(pheno$study_id,fam1$V2)]-1
pheno$female[which(!pheno$female %in% c(0,1))]=NA
pheno$case_status = 1

pheno_full = pheno[complete.cases(pheno),]

write.table(pheno_full,file='intdata/08-BIOJUME_JME_ids.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
```

In bash, reformat data:

```{bash}
# Subset BIOJUME data
plink2 --bfile rawdata/BIOJUME_R1234567 --keep intdata/08-BIOJUME_JME_ids.txt --make-pgen --out intdata/08-BIOJUME_JME

# Replace hyphens with underscores in psam file
sed 's/-/_/g' intdata/08-BIOJUME_JME.psam > tmp.txt
mv tmp.txt intdata/08-BIOJUME_JME.psam

# Prepare and add phenotypes
awk -F'\t' '{print $1 "\t" $2 "\t2"}' intdata/08-BIOJUME_JME.psam > intdata/08-BIOJUME_phenotypes.txt

awk -F'\t' 'NR==1 {$3="STATUS"} {print $0}' OFS='\t' intdata/08-BIOJUME_phenotypes.txt > temp.txt
mv temp.txt intdata/08-BIOJUME_phenotypes.txt

# Create file ready to merge with SFS data for PCA:
plink2 --pfile intdata/08-BIOJUME_JME --pheno intdata/08-BIOJUME_phenotypes.txt --make-pgen --out intdata/09-BIOJUME_prepared_genotype

# Final tweaks to allow merge

awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $5, $6}' intdata/09-BIOJUME_prepared_genotype.psam > tidied.psam
mv tidied.psam intdata/09-BIOJUME_prepared_genotype.psam

cut -f 1-5 intdata/09-BIOJUME_prepared_genotype.pvar > tidied.pvar
mv tidied.pvar intdata/09-BIOJUME_prepared_genotype.pvar
```

## Subset and Format Imputed Data

In bash, merge VCF files and convert to plink format.

```{bash}
# Make sex file that will work using vcf-ID-format
awk 'NR > 1 { gsub(/_/, "-", $2); gsub(/_/, "-", $5); print $2 "\t" $5 }' intdata/09-BIOJUME_prepared_genotype.psam > intdata/08-BIOJUME_imputed_sexupdate.txt

# Convert each chromosome file to plink format
for chr in {1..22} X; do
  echo "Beginning Chr $chr"
  
  # Specify the path to original VCF file
  vcf_file="rawdata/BIOJUME_R1234567.rsq08.chr${chr}.vcf.gz"
  
  # Specify the path for the updated VCF file
  output_file="intdata/08-BIOJUME_imputed_chr${chr}.dose.vcf"
  
  # Extract the header line from the original VCF file
  header=$(zgrep "^#" "$vcf_file")
    
  echo "Replacing..."
  
  # Replace underscores with hyphens in the sample IDs of the header line
  updated_header=$(echo "$header" | awk -F'\t' -v OFS='\t' '{
    if ($0 ~ /^#CHROM/) {
        for (i = 10; i <= NF; i++) {
            sub(/Round[^_]+_/, "&-", $i)
           gsub("_", "-", $i)
            gsub("--", "_", $i)
        }
   }
    print
  }')
  
  echo "Saving output..."
  # Output the updated header line
  echo "$updated_header" > "$output_file"

  # Output the rest of the VCF file (excluding the original header)
  zgrep -v "^#" "$vcf_file" >> "$output_file"

  plink2 --vcf intdata/08-BIOJUME_imputed_chr"${chr}".dose.vcf dosage=DS --update-sex intdata/08-BIOJUME_imputed_sexupdate.txt --make-pgen --out intdata/08-BIOJUME_imputed_chr"${chr}"
done

# Merge plink files
for chr in {1..22} X; do
  echo "intdata/08-BIOJUME_imputed_chr${chr}" >> intdata/08-BIOJUME_imputed_merge_list.txt
done

# before merge have to replace duplicated rsids in each file with CHR:POS:REF:ALT ID names
for chr in {1..22} X; do
  # Find duplicate IDs
  awk '!/^##/ && !/^#CHROM/ {print $3}' intdata/08-BIOJUME_imputed_chr${chr}.pvar | sort | uniq -d > dup_ids.txt
  
  # ifffff there are duplicate ID's, rename with CHR:POS:REF:ALT
  awk -v dups=dup_ids.txt '
BEGIN {
  while (getline line < dups) dup[line]=1
}
# Skip metadata and header
/^##/ || /^#CHROM/ {print; next}
{
  if ($3 in dup) {       # duplicate rsID
    $3 = $1":"$2":"$4":"$5  # CHR:POS:REF:ALT
  }
  print
}' intdata/08-BIOJUME_imputed_chr${chr}.pvar > intdata/08-BIOJUME_imputed_chr${chr}_fixed.pvar

  mv intdata/08-BIOJUME_imputed_chr${chr}_fixed.pvar intdata/08-BIOJUME_imputed_chr${chr}.pvar
done

plink2 --pmerge-list intdata/08-BIOJUME_imputed_merge_list.txt --make-pgen --out intdata/08-BIOJUME_all_imputed

rm intdata/08-BIOJUME_imputed*
```

In R, adjust the formatting of the psam file.

```{r}
psam = read.table('intdata/08-BIOJUME_all_imputed.psam',header=FALSE)

psam$V1 = gsub('-','_',psam$V1)

psam2 = data.frame(V1=psam$V1,"IID"=psam$V1,"SEX"=psam$V2)
colnames(psam2)[1]="#FID"

write.table(psam2,file=paste0('intdata/08-BIOJUME_all_imputed.psam'),row.names=FALSE,col.names=TRUE,quote=FALSE)
```

In bash, subset to only the probable/definite JME cases and adjust format of psam/pvar files further:

```{bash}
plink2 \
  --pfile intdata/08-BIOJUME_all_imputed \
  --keep intdata/09-BIOJUME_prepared_genotype.psam \
  --make-pgen \
  --out intdata/09-BIOJUME_prepared_imputed
  
# Use the genotype psam to make set of FID_IID keys, use those keys to get PAT/MAT/SEX/STATUS values for the imputed data psam file
awk '
NR==FNR {
    if (NR>1) { key=$1"_"$2; line[key]=$0 }
    next
}

FNR==1 {
    # Replace header with full header from set 2
    print "#FID IID PAT MAT SEX STATUS"
    next
}

{
    key=$1"_"$2
    if (key in line) {
        print line[key]
    } else {
        # If a sample is missing in reference, keep the original minimal info and fill missing fields
        print $1, $2, "0", "0", "0", "0"
    }
}
' intdata/09-BIOJUME_prepared_genotype.psam intdata/09-BIOJUME_prepared_imputed.psam \
> intdata/09-BIOJUME_prepared_imputed_updated.psam

# Check that row order was mantained
tail intdata/09-BIOJUME_prepared_imputed_updated.psam
tail intdata/09-BIOJUME_prepared_imputed.psam

mv intdata/09-BIOJUME_prepared_imputed_updated.psam intdata/09-BIOJUME_prepared_imputed.psam

awk 'BEGIN{FS=OFS="\t"} 
NR==1 {print "#FID", "IID", "SEX", "STATUS"} 
NR>1 {print $1, $2, $5, $6}' intdata/09-BIOJUME_prepared_imputed.psam > tidied.psam
mv tidied.psam intdata/09-BIOJUME_prepared_imputed.psam

tail -n +50 intdata/09-BIOJUME_prepared_imputed.pvar | cut -f 1-5 > tidied.pvar
mv tidied.pvar intdata/09-BIOJUME_prepared_imputed.pvar
```

