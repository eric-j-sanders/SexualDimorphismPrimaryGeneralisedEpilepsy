# Purpose

When the codes in this folder are executed, three secondary analytical steps are completed to prepare additional SNP filters for application to the GWAS results:

* HWE testing specifically in controls (SNPs will be removed if they have P<1e-4 in this test)
* A control vs. control GWAS where both groups are of European ancestry (SNPs will be removed if they have P<1e-5 in this GWAS)
* Chi-square testing of AF differences between completing imputation separately by batch or in merged batches (SNPs will be removed if both imputation pipelines have imputation R2>0.8 but the SNP has P<1e-4 in the chi-square test of AF difference)

# Input files

* The listed outputs of all previous repository blocks
* A set of psam/pgen/pvar files in plink format, documenting CF sample genotype data (before imputation), saved to relative path rawdata/05-all_CF_genotype_clean. This is after QC steps previously detailed.
* A set of psam/pgen/pvar files in plink format, documenting CF sample imputed genotype data, saved to relative path rawdata/05-all_CF_imputed_clean. This is after QC steps previously detailed.
* The scripts saved to this directory
* An "extended_ld_regions.txt" file that lists extended LD regions (HG38) for removal before completing PCA.

# Primary outputs

* A list of SNPs that pass the HWE test in controls (P > 1e-4) at relative path intdata/29-control_HWE_pass_e04_snps.txt
* A set of GWAS summary statistics satisfying P<0.1 from a SFS vs. CF GWAS, saved in the relative path intdata/47-gwas_out_noint_G_pfilt.txt
* A list of SNPs that had significant AF difference when imputed by batch or jointly, despite having imputation R2>0.8 in both pipelines, is saved to the relative path intdata/86-variants_with_differing_af_by_pipeline_CASES_OR_CONTROLS.txt

# Code

## HWE Testing in Controls

Via bash:

```{bash}
plink2 --pfile intdata/14-merge_imputed_eur_unrel --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --hardy --out intdata/29-control_HWE

awk -F'\t' 'NR>1 && $10>1e-4 {print $2}' intdata/29-control_HWE.hardy > intdata/29-control_HWE_pass_e04_snps.txt

awk -F'\t' 'NR>1 && $14>1e-4 {print $2}' intdata/29-control_HWE.hardy.x >> intdata/29-control_HWE_pass_e04_snps.txt
```

## Secondary GWAS Between SFS and CF Controls

### Ancestry in CF Sample

Ancestry analysis of CF controls:

```{bash}
plink2 --pfile rawdata/05-all_CF_genotype_clean --make-bed --out intdata/36-CF_forgrafpop

grafpop intdata/36-CF_forgrafpop.bed intdata/36-CF-grafpop.txt

perl ~/Programs/PlotGrafPopResults.pl intdata/36-CF-grafpop.txt figures/CF_GRAFPOP_output.png

perl ~/Programs/SaveSamples.pl intdata/36-CF-grafpop.txt intdata/36-CF-grafpop_samples.txt

rm intdata/36-CF_forgrafpop*
```

Using R, a more tidy file is created that marks FID, IID, inferred population for a combined set of CF and SFS control samples:

```{R}
psam = read.table('rawdata/05-all_CF_genotype_clean.psam')

grafpop = read.table('intdata/36-CF-grafpop_samples.txt',comment.char='',header=TRUE,sep='\t')

merged = merge(psam,grafpop,by.x='V2',by.y='Subject')

merged$inferred_european = as.numeric(merged$Computed.population=='European')

merged = merged[,c('V1','V2','PopID','Computed.population','inferred_european')]

colnames(merged)[1:2]=c('FID',"IID")

write.csv(merged,file='intdata/36-CF_ancestry.csv',col.names=TRUE,row.names=FALSE)

sfs = read.csv('intdata/10-SFS2_ancestry.csv')
cf = read.csv('intdata/36-CF_ancestry.csv')

eur_ids = rbind(sfs[which(sfs$inferred_european==1),c(1,2)],cf[which(cf$inferred_european==1),c(1,2)])

write.table(eur_ids,file='intdata/37-eur_cf_sfs_ids.txt',sep='\t',row.names=FALSE,quote=FALSE)
```

### Merging CF and SFS, subsetting to Europeans

Genotype datasets are merged in bash:

```{bash}
sbatch code/merge_bplink_files_v3_for_plink2_whenFile2HasChrPosRefAlt.sh intdata/07-SFS2_prepared_genotype_nopsych rawdata/05-all_CF_genotype_clean intdata/38-merge_cf_sfs_geno genotypemerge
```

Using R, genotype psam is cleaned up post-merge:

```{R}
psam = read.table('intdata/38-merge_cf_sfs_geno.psam',comment.char='',header=TRUE)
psam$`#FID` = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][1])
psam$IID = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][2])

psam$`#FID` = gsub('-','_',psam$`#FID`)
psam$IID = gsub('-','_',psam$IID)

psam2 = psam[,c('#FID','IID','SEX')]

psam_cf = read.table('rawdata/05-all_CF_genotype_clean.psam',comment.char='',header=TRUE)
colnames(psam_cf)[1]='#FID'
psam_cf$STATUS=1

psam_sfs = read.table('intdata/07-SFS2_prepared_genotype_nopsych.psam',comment.char='',header=TRUE)
colnames(psam_sfs)[1]='#FID'
psam_sfs$STATUS=2

psam_ogs = rbind(psam_cf,psam_sfs)
psam_ogs$IID = gsub('-','_',psam_ogs$IID)
psam_ogs$`#FID` = gsub('-','_',psam_ogs$`#FID`)

psam2$STATUS=NA

psam2$SEX[match(psam_ogs$IID, psam2$IID)] = psam_ogs$SEX
psam2$STATUS[match(psam_ogs$IID, psam2$IID)] = psam_ogs$STATUS

write.table(psam2,'intdata/38-merge_cf_sfs_geno.psam',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
```

Using bash, the genotype data is subset to Europeans:

```{bash}
sed 's/-/_/g' intdata/37-eur_cf_sfs_ids.txt > fix.txt
mv fix.txt intdata/37-eur_cf_sfs_ids.txt

plink2 --pfile intdata/38-merge_cf_sfs_geno --keep intdata/37-eur_cf_sfs_ids.txt --make-pgen --out intdata/39-merge_cf_sfs_genotype_eur
```

And the imputed data is also merged and subset:

```{bash}
sbatch code/merge_bplink_files_v3_for_plink2_whenFile2HasChrPosRefAlt.sh intdata/07-SFS2_prepared_imputed_nopsych rawdata/05-all_CF_imputed_clean intdata/38-merge_cf_sfs_imp genmerge
```

And using R, the psam is fixed post-merge:

```{r}
psam = read.table('intdata/38-merge_cf_sfs_imp.psam',comment.char='',header=TRUE)

psam$`#FID` = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][1])
psam$IID = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][2])

psam$`#FID` = gsub('-','_',psam$`#FID`)
psam$IID = gsub('-','_',psam$IID)

psam2 = psam[,c('#FID','IID','SEX')]

psam_cf = read.table('rawdata/05-all_CF_imputed_clean.psam',comment.char='',header=TRUE)
colnames(psam_cf)[1]='#FID'
psam_cf$STATUS=2
psam_cf$`#FID` = gsub('-','_',psam_cf$`#FID`)
psam_cf$IID = gsub('-','_',psam_cf$IID)

psam_sfs = read.table('intdata/07-SFS2_prepared_imputed_nopsych.psam',comment.char='',header=TRUE)
colnames(psam_sfs)[1]='#FID'
psam_sfs$`#FID` = gsub('-','_',psam_sfs$`#FID`)
psam_sfs$IID = gsub('-','_',psam_sfs$IID)

psam_ogs = rbind(psam_cf,psam_sfs)

psam2$STATUS=NA

psam2$SEX[match(psam_ogs$IID, psam2$IID)] = psam_ogs$SEX
psam2$STATUS[match(psam_ogs$IID, psam2$IID)] = psam_ogs$STATUS

write.table(psam2,'intdata/38-merge_cf_sfs_imp.psam',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
```

And finally subset to Europeans:

```{bash}
sbatch --mem=120G --cpus-per-task=4 --time=04:00:00 --
wrap="plink2 --pfile intdata/38-merge_cf_sfs_imp --keep intdata/37-eur_cf_sfs_ids.txt --make-pgen --out intdata/39-merge_cf_
sfs_imp_eur"
```

### Relatedness and PCA in SFS-CF population

In bash:

```{bash}
plink2 --pfile intdata/39-merge_cf_sfs_genotype_eur --make-bed --out intdata/39-forking
plink --bfile intdata/39-forking --make-bed --out intdata/39-forking2

king -b intdata/39-forking2.bed --related --degree 2 --prefix intdata/39-kinship_all
king -b intdata/39-forking2.bed --unrelated --degree 2 --prefix intdata/39-kinship_unrel

rm intdata/39-forking*

plink2 --pfile intdata/39-merge_cf_sfs_genotype_eur --keep intdata/39-kinship_unrelunrelated.txt --make-pgen --out intdata/40-merge_genotype_eur_unrel
```

Use result to subset the imputed data to unrelated samples:

```{bash}
sbatch --mem=120G --cpus-per-task=4 --time=04:00:00 --wrap="plink2 --pfile intdata/39-merge_cf_sfs_imp_eur --keep intdata/39-kinship_unrelunrelated.txt --make-pgen --out intdata/40-merge_imp_eur_unrel"
```

Removing Chr7 before continuing:

```{bash}
sbatch --mem=120G --cpus-per-task=4 --time=04:00:00 --wrap="plink2 --pfile intdata/40-merge_imp_eur_unrel --not-chr 7 --make-pgen --out intdata/40-merge_imp_eur_unrel_noCFTR"
```

Performing PCA on the genotype data:

```{bash}
bash code/prePCA_plink2.sh intdata/40-merge_genotype_eur_unrel_noCFTR intdata/41

bash code/plink2_PCA.sh 'intdata/41-pruned' 'intdata/42-pca_res'
```

Using R a table of PCA results is assembled:

```{r}
psam = read.table('intdata/40-merge_genotype_eur_unrel_noCFTR.psam')
colnames(psam)=c('FID','IID','SEX','STATUS')

eigenval = read.table('intdata/42-pca_res.eigenval')

eigenvec = read.table('intdata/42-pca_res.eigenvec')
colnames(eigenvec) = c('FID','IID',paste0('PC',1:(ncol(eigenvec)-2)))

outfile = 'intdata/43-eur_unrel_pheno.txt'

eigenval$V2 = 1:nrow(eigenval)
eigenval$V3 = eigenval$V1/sum(eigenval$V1)

# Save first 10 PCs with SEX+PHENO in new file
new_pheno = merge(psam, eigenvec, by = "IID", all = FALSE)[,c(2,1,3,4,6:15)]
colnames(new_pheno)[1] = 'FID'

new_pheno$Pheno = factor(new_pheno$STATUS,levels=c(1,2),labels=c('CF','SFS'))

new_pheno$IID = gsub('-','_',new_pheno$IID)
new_pheno$FID = gsub('-','_',new_pheno$FID)

write.table(new_pheno[,1:14], outfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
```

### Case-Control Matching

In R:

```{r}
library(PCAmatchR)
library(optmatch)

options("optmatch_max_problem_size" = Inf)

cov = read.table('intdata/43-eur_unrel_pheno.txt',header=TRUE)

eigen = read.table("intdata/42-pca_res.eigenval")$V1

cov$formatch = as.numeric(cov$STATUS==2)

k=10

matches = match_maker(PC=cov[,c(2,5:(k+4))],eigen_value=eigen[1:k]/sum(eigen[1:k]),data=cov,ids='IID',case_control='formatch',num_controls=1,num_PCs=k)

write.table(matches$matches[,c(2,1,5:(k+4))],file=paste0('intdata/44-cov_matched_',k,'.txt'),row.names=FALSE,col.names=TRUE,quote=FALSE)
```

Subset via bash to the matched samples with hardcall genotype format:

```{bash}
sbatch --mem=120G --cpus-per-task=4 --time=04:00:00 --wrap="plink2 --pfile intdata/40-merge_imp_eur_unrel_noCFTR --keep intdata/44-cov_matched_10.txt --hard-call-threshold 0.49999 --make-bed --out intdata/45-merge_imputed_eur_unrel_matched"
```

### Running CF-SFS GWAS

Final covariate file assembled in R:

```{r}
cov = read.table('intdata/44-cov_matched_10.txt',header=TRUE)
psam = read.table('intdata/45-merge_imputed_eur_unrel_matched.psam',header=TRUE,comment='')

covar = merge(cov,psam[,c('IID','SEX')],by='IID')
covar = covar[,c(2,1,3:ncol(covar))]

write.table(covar,file='intdata/45-covar_for_GWAS.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

colnames(psam)=c('#FID','IID','SEX',"STATUS")

write.table(psam[,c(1,2,4)],'intdata/45-merge_imputed_eur_unrel_matched.psam',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
```

And the GWAS is completed via bash:

```{bash}
sbatch --mem=120G --cpus-per-task=4 --time=03:00:00 --wrap="plink2 \
  --bfile intdata/45-merge_imputed_eur_unrel_matched \
  --maf 0.01 \
  --glm no-x-sex\
  --covar intdata/45-covar_for_GWAS.txt \
  --out intdata/46-GWAS_res_noint"
```

In bash the GWAS results are processed into a more tidy format (the code was adapted from previous code that processed several GWAS results sets -- this is why it "loops" through a single set of parameters)

```{bash}
input_files=("intdata/46-GWAS_res_noint.PHENO1.glm.logistic.hybrid")

main_test=("ADD")

mid_files=("intdata/47-gwas_out_noint_G.txt")

sig_files=("intdata/47-gwas_out_noint_G_pfilt.txt")

for i in "${!input_files[@]}"; do

  input_file="${input_files[i]}"
  mid_file="${mid_files[i]}"
  sig_file="${sig_files[i]}"
  test="${main_test[i]}"

  echo "Looking for ${test} in ${input_file}"

  temp_file=$(mktemp)

  awk '$16 != "NA"' $input_file > $mid_file

  head -n 1 $mid_file > $temp_file
  awk '$11 == test' test="$test" "$mid_file" >> "$temp_file"
  
  mv "$temp_file" "$mid_file"
  
  temp_file=$(mktemp)

  awk '$16 != "NA"' $input_file > $sig_file

  head -n 1 $sig_file > $temp_file
  awk '$16 < 0.1 && $11 == test' test="$test" "$sig_file" >> "$temp_file"
  
  mv "$temp_file" "$sig_file"

  # These files don't require the columns 6,8,10,12,14,15,17:

  temp_file=$(mktemp)

  # Use cut to remove columns 6,8,10,12,14,15,17 and save the result to the temporary file
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5 ,$7, $9, $11, $13, $16}' "$mid_file" > "$temp_file"

  # Move the temporary file to overwrite the original file
  mv "$temp_file" "$mid_file"

  temp_file=$(mktemp)

  # Use cut to remove columns 6,8,10,12,14,15,17 and save the result to the temporary file
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5 ,$7, $9, $11, $13, $16}' "$sig_file" > "$temp_file"

  # Move the temporary file to overwrite the original file
  mv "$temp_file" "$sig_file"
  
done
```

## Comparison of Imputation Pipeline (Identifying SNPs with significant AF variation)

### Re-running genotype imputation on a merged sample

To contrast the main GWAS pipeline in which individual sample batches were already imputed during previously described QC steps, in this section the imputation was repeated, with the only difference being that the samples were merged across batches first.

First, the unrelated European merged set of cases and controls are formatted into VCFs that can be uploaded to TopMed, using bash:

```{bash}
for chr in {1..22}; do
  echo "Beginning Chr $chr"

  # Export VCF from plink
  plink2 --pfile intdata/14-merge_genotype_eur_unrel --chr $chr --recode vcf --out intdata/70-chr${chr}

  # Rename chromosomes, sort, compress
  bcftools annotate --rename-chrs <(awk 'BEGIN{for(i=1;i<=22;i++) print i,"\tchr"i; print "X\tchrX"}') intdata/70-chr${chr}.vcf | bcftools sort -Oz -o intdata/71-chr${chr}_topmed.vcf.gz

  # Change header version
  bcftools view -h intdata/71-chr${chr}_topmed.vcf.gz | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > header.txt

  bcftools reheader -h header.txt intdata/71-chr${chr}_topmed.vcf.gz -o intdata/71-chr${chr}_topmed.vcf.gz.tmp

  mv intdata/71-chr${chr}_topmed.vcf.gz.tmp intdata/71-chr${chr}_topmed.vcf.gz

  # Index
  bcftools index -t intdata/71-chr${chr}_topmed.vcf.gz

  rm intdata/70-chr${chr}.vcf header.txt

done

#ChrX is more involved, must split and fix het hap calls

chr='X'

echo "Beginning Chr $chr"

plink2 --pfile intdata/14-merge_genotype_eur_unrel --set-hh-missing --make-pgen --out temp_clean

# Export VCF from plink
plink2 --pfile temp_clean --chr $chr --recode vcf --out intdata/70-chr${chr}

# Rename chromosomes, sort, compress
bcftools annotate --rename-chrs <(awk 'BEGIN{for(i=1;i<=22;i++) print i,"\tchr"i; print "X\tchrX"}') intdata/70-chr${chr}.vcf | bcftools sort -Oz -o intdata/71-chr${chr}_topmed.vcf.gz

# Change header version
bcftools view -h intdata/71-chr${chr}_topmed.vcf.gz | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > header.txt

bcftools reheader -h header.txt intdata/71-chr${chr}_topmed.vcf.gz -o intdata/71-chr${chr}_topmed.vcf.gz.tmp

mv intdata/71-chr${chr}_topmed.vcf.gz.tmp intdata/71-chr${chr}_topmed.vcf.gz

# Index
bcftools index -t intdata/71-chr${chr}_topmed.vcf.gz

rm intdata/70-chr${chr}.vcf header.txt
rm temp*
```

At this point the VCFs were uploaded to TopMed imputation server, using default parameters and no R2 filter.

The imputation results are saved to the folder intdata/72-imputation_res and unzipped using the provided password (code not shown)

Using bash, the imputed data are converted to plink format, and subset to R2>0.8:

```{bash}
for chr in {1..22}; do
  echo "Beginning Chr $chr"
  
  # Specify the path to original VCF file
  vcf_file="intdata/72-imputation_res/chr${chr}.dose.vcf.gz"
  
  # Fix sample IDs and filter down to R2 thresholds
  
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=filter_chr${chr}
#SBATCH --output=job_output/filter_chr${chr}.log
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00

export PATH=$PATH:$SLURM_SUBMIT_DIR
source ~/Programs/miniconda3/etc/profile.d/conda.sh
conda activate gwas2024

threads=\$SLURM_CPUS_PER_TASK

# Filter for R2 >= 0.3
bcftools view --threads \$threads -i 'R2>=0.3' "$vcf_file" -Oz -o "temp_${chr}_${SLURM_JOB_ID}.vcf.gz" &&\
plink2 --threads \$threads --vcf "temp_${chr}_${SLURM_JOB_ID}.vcf.gz" dosage=DS --double-id --make-pgen --out intdata/73-joint_imputed_r3_chr"${chr}" &&\

# Filter the temp file for R2 >= 0.8
bcftools view --threads \$threads -i 'R2>=0.8' "temp_${chr}_${SLURM_JOB_ID}.vcf.gz" -Oz -o "temp2_${chr}_${SLURM_JOB_ID}.vcf.gz" &&\
plink2 --threads \$threads --vcf "temp2_${chr}_${SLURM_JOB_ID}.vcf.gz" dosage=DS --double-id --make-pgen --out intdata/73-joint_imputed_r8_chr"${chr}" &&\

# Clean up
rm "temp_${chr}_${SLURM_JOB_ID}.vcf.gz" "temp2_${chr}_${SLURM_JOB_ID}.vcf.gz"
EOF

done
```

As a rough comparison: 

* When imputed separately: the SFS batch took 0.5m genotyped SNPs and produced 10.9m imputed with R2>0.8, and individual BIOJUME rounds took 0.5m-1.8m SNPs in and returned 10m-22m SNPs per round with R2>0.8, although the final intersecting set was 7.3m with R2>0.8 in all batches.
* Now after doing this joint imputation, 157,222 Total SNPs were uploaded and 25 million were imputed with R2>0.8

Post-imputation data cleanup in R:

```{r}
library(stringr)
library(dplyr)

cov=read.table('intdata/17-eur_unrel_pheno.txt',header=TRUE)

for(i in 1:22){ 
    for(r in c(3,8)){
        fam = read.table(paste0('intdata/73-joint_imputed_r',r,'_chr',i,'.psam'),header=TRUE,comment.char='')
      colnames(fam)=c('#FID','IID','SEX')
      fam=fam[,c(1,2)]
      counts = str_count(fam$IID,'_')
      
      fam$IID[which(counts==1)]=vapply(fam$IID[which(counts==1)],function(x) strsplit(x,'_')[[1]][2], FUN.VALUE=character(1))
      fam$`#FID`[which(counts==1)]=vapply(fam$`#FID`[which(counts==1)],function(x) strsplit(x,'_')[[1]][1], FUN.VALUE=character(1))
      
      fam$IID[which(counts==3)]=vapply(fam$IID[which(counts==3)],function(x) paste0(strsplit(x,'_')[[1]][c(3,4)],collapse='_'), FUN.VALUE=character(1))
      fam$`#FID`[which(counts==3)]=vapply(fam$`#FID`[which(counts==3)],function(x) paste0(strsplit(x,'_')[[1]][c(1,2)],collapse='_'), FUN.VALUE=character(1))
      
      fam = left_join(fam,cov[,c(1:4)],by=c('#FID'='FID','IID'='IID'))
      
      write.table(fam,file=paste0('intdata/73-joint_imputed_r',r,'_chr',i,'.psam'),row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
  
}
```

### Preparing New Imputed Dataset for Comparison with Original GWAS Dataset

Using bash, apply MAF>0.01 filter to both cases and controls.

```{bash}
for chr in {1..22}; do
    plink2 --pfile "intdata/73-joint_imputed_r8_chr${chr}" --keep intdata/09-BIOJUME_prepared_imputed.psam --freq counts --out "intdata/74-cases_freq_${chr}"
    
    plink2 --pfile "intdata/73-joint_imputed_r8_chr${chr}" --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --freq counts --out "intdata/74-controls_freq_${chr}"
    
    awk 'NR>1 {
      alt=$5; obs=$6;
      maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
      if (maf > 0.01) print $2
    }' "intdata/74-controls_freq_${chr}.acount" > "intdata/75-control_snps_${chr}_maf01.txt"
    
    awk 'NR>1 {
      alt=$5; obs=$6;
      maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
      if (maf > 0.01) print $2
    }' "intdata/74-cases_freq_${chr}.acount" > "intdata/75-case_snps_${chr}_maf01.txt"
    
    # Only keep SNPs with MAF>0.01 in both cases and controls
    comm -12 <(sort "intdata/75-case_snps_${chr}_maf01.txt") <(sort "intdata/75-control_snps_${chr}_maf01.txt") > "intdata/75-snps_maf01_${chr}_overlap.txt"
    
    sbatch --mem=80G --cpus-per-task=2 --time=03:00:00 --wrap="plink2 --pfile \"intdata/73-joint_imputed_r8_chr${chr}\" --extract \"intdata/75-snps_maf01_${chr}_overlap.txt\" --make-pgen --out \"intdata/76-joint_imputed_r8_chr${chr}_maf01both\""
done
```

Merge chromosomes:

```{bash}
for chr in {1..22}; do
  echo "intdata/76-joint_imputed_r8_chr${chr}_maf01both" >> intdata/77-merge_list.txt
done

sbatch --mem=100G --cpus-per-task=2 --time=06:00:00 --job-name="PLNK_MERGE" --wrap='plink2 --pmerge-list intdata/77-merge_list.txt --max-alleles 2 --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 55 --make-pgen --out intdata/78-final_merged_imputed'
```

Do control-specific HWE control and round to hard calls:

```{bash}
plink2 --pfile intdata/78-final_merged_imputed --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --hardy --out intdata/79-control_HWE

awk -F'\t' 'NR>1 && $10>1e-4 {print $2}' intdata/79-control_HWE.hardy > intdata/79-control_HWE_pass_e04_snps.txt

plink2 --pfile intdata/78-final_merged_imputed --extract intdata/79-control_HWE_pass_e04_snps.txt --hard-call-threshold 0.49999 --make-bed --out intdata/80-joint_imputed_maf01both_hwe04_hardcalls
```

### Compare AF between Imputation Pipelines

Separately for cases and controls, collect AF information for SNPs from both the original by-batch imputation pipeline and this new pre-merged imputation pipeline.

```{bash}
# First in cases

plink2 --pfile intdata/14-merge_imputed_eur_unrel --keep intdata/09-BIOJUME_prepared_imputed.psam --freq cols=+pos --out intdata/86-freq_batchimpute_cases

plink2 --bfile intdata/80-joint_imputed_maf01both_hwe04_hardcalls --freq cols=+pos --keep intdata/09-BIOJUME_prepared_imputed.psam --out intdata/86-freq_jointimpute_cases

# Second in controls
plink2 --pfile intdata/14-merge_imputed_eur_unrel --freq cols=+pos --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --out intdata/86-freq_batchimpute_controls

plink2 --bfile intdata/80-joint_imputed_maf01both_hwe04_hardcalls --freq cols=+pos --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --out intdata/86-freq_jointimpute_controls
```

And in R complete chi-square tests of differences in AF, for either cases or for controls:

```{r}
sbatch --mem=30G --cpus-per-task=2 --time=04:00:00 --wrap="Rscript code/pipeline_chisq.R"
```

This creates the output file intdata/86-variants_with_differing_af_by_pipeline_CASES_OR_CONTROLS.txt that marks SNPs to remove because of imputation pipeline being highly influential to imputed AF, despite both pipelines appearing "confident"