# Purpose

When the codes in this folder are executed, a sex-stratified GWAS of JME is completed using logistic regression in plink. Regression results are processed to filter out SNPs based on previously implemented secondary steps.

# Input files

* The listed outputs of all previous repository blocks

# Primary outputs

* Pre-GWAS genetic datasets in bed/bim/fam plink format for each sex, saved with prefixes intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls and intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls
* GWAS summary statistics from the female-only case control analysis, saved to intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt
* GWAS summary statistics from the male-only case control analysis, saved to intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt
* Mirrored Manhattan plot comparing the two GWAS results: figures/my_miami_plot_refined.png

# Code

## Sex-Specific Case-Control Matching

Using R, cases are matched to controls, first in males and then in females.

```{r}
library(PCAmatchR)
library(optmatch)

cov = read.table('intdata/17-eur_unrel_pheno.txt',header=TRUE)

eigen = read.table("intdata/16-pca_res.eigenval")$V1

cov$formatch = as.numeric(cov$STATUS==2)

cov$SEX = as.factor(cov$SEX)

k=10

matches_m = match_maker(PC=cov[cov$SEX==1,c(2,5:(k+4))],eigen_value=eigen[1:k]/sum(eigen[1:k]),data=cov[cov$SEX==1,],ids='IID',case_control='formatch',num_controls=5,num_PCs=k)

matches_f = match_maker(PC=cov[cov$SEX==2,c(2,5:(k+4))],eigen_value=eigen[1:k]/sum(eigen[1:k]),data=cov[cov$SEX==2,],ids='IID',case_control='formatch',num_controls=2,num_PCs=k)

newset_m = matches_m$matches[,c(2,1,5:(k+4))]
newset_f = matches_f$matches[,c(2,1,5:(k+4))]

write.table(newset_m,file=paste0('intdata/51-cov_matched_males_',k,'.txt'),row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(newset_f,file=paste0('intdata/51-cov_matched_females_',k,'.txt'),row.names=FALSE,col.names=TRUE,quote=FALSE)
```

Male-specific GWAS will be 236 cases vs. 1180 controls. (5:1 matching)

Female-specific GWAS will be 473 cases vs. 946 controls. (2:1 matching)

## Completion of Initial GWAS

Prepare GWAS sets by subsetting to match sets, then applying MAF>0.01 filter to both cases and controls:

```{bash}
plink2 --pfile intdata/14-merge_imputed_eur_unrel --keep intdata/51-cov_matched_males_10.txt --make-pgen --out intdata/52-merge_males_imputed_eur_unrel_matched

plink2 --pfile intdata/14-merge_imputed_eur_unrel --keep intdata/51-cov_matched_females_10.txt --make-pgen --out intdata/52-merge_females_imputed_eur_unrel_matched

#################################
# Reduce male SNPs by MAF first:

plink2 --pfile intdata/52-merge_males_imputed_eur_unrel_matched --keep intdata/09-BIOJUME_prepared_imputed.psam --freq counts --out intdata/52-males_cases_freq

plink2 --pfile intdata/52-merge_males_imputed_eur_unrel_matched --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --freq counts --out intdata/52-males_controls_freq

awk 'NR>1 {
  alt=$5; obs=$6;
  maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
  if (maf > 0.01) print $2
}' intdata/52-males_cases_freq.acount > intdata/52-males_control_snps_maf01.txt

awk 'NR>1 {
  alt=$5; obs=$6;
  maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
  if (maf > 0.01) print $2
}' intdata/52-males_controls_freq.acount > intdata/52-males_case_snps_maf01.txt

# Only keep SNPs with MAF>0.01 in both cases and controls
comm -12 <(sort intdata/52-males_case_snps_maf01.txt) <(sort intdata/52-males_control_snps_maf01.txt) > intdata/53-males_snps_maf01_overlap.txt

plink2 --pfile intdata/52-merge_males_imputed_eur_unrel_matched --extract intdata/53-males_snps_maf01_overlap.txt --make-pgen --out intdata/54-merge_males_imputed_eur_unrel_matched_maf01both

#################################
# Now reduce female SNPs by MAF:

plink2 --pfile intdata/52-merge_females_imputed_eur_unrel_matched --keep intdata/09-BIOJUME_prepared_imputed.psam --freq counts --out intdata/52-females_cases_freq

plink2 --pfile intdata/52-merge_females_imputed_eur_unrel_matched --keep intdata/07-SFS2_prepared_imputed_nopsych.psam --freq counts --out intdata/52-females_controls_freq

awk 'NR>1 {
  alt=$5; obs=$6;
  maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
  if (maf > 0.01) print $2
}' intdata/52-females_cases_freq.acount > intdata/52-females_control_snps_maf01.txt

awk 'NR>1 {
  alt=$5; obs=$6;
  maf = (alt <= obs-alt ? alt/obs : (obs-alt)/obs);
  if (maf > 0.01) print $2
}' intdata/52-females_controls_freq.acount > intdata/52-females_case_snps_maf01.txt

# Only keep SNPs with MAF>0.01 in both cases and controls
comm -12 <(sort intdata/52-females_case_snps_maf01.txt) <(sort intdata/52-females_control_snps_maf01.txt) > intdata/53-females_snps_maf01_overlap.txt

plink2 --pfile intdata/52-merge_females_imputed_eur_unrel_matched --extract intdata/53-females_snps_maf01_overlap.txt --make-pgen --out intdata/54-merge_females_imputed_eur_unrel_matched_maf01both
```

Round genotypes to hard 0/1/2 calls and subset to SNPs with HWE>1e-4 in controls:

```{bash}
plink2 --pfile intdata/54-merge_females_imputed_eur_unrel_matched_maf01both --extract intdata/29-control_HWE_pass_e04_snps.txt --hard-call-threshold 0.49999 --make-bed --out intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls

plink2 --pfile intdata/54-merge_males_imputed_eur_unrel_matched_maf01both --extract intdata/29-control_HWE_pass_e04_snps.txt --hard-call-threshold 0.49999 --make-bed --out intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls
```

And run the two GWAS:

```{bash}
sbatch --mem=80G --cpus-per-task=2 --time=03:00:00 --wrap="plink2 \
  --bfile intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
  --glm \
  --covar intdata/51-cov_matched_males_10.txt \
  --out intdata/56-GWAS_res_males_hardcall"
  
sbatch --mem=80G --cpus-per-task=2 --time=03:00:00 --wrap="plink2 \
  --bfile intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
  --glm \
  --covar intdata/51-cov_matched_females_10.txt \
  --out intdata/56-GWAS_res_females_hardcall"
```

And tidy the results:

```{bash}
input_files=("intdata/56-GWAS_res_females_hardcall.PHENO1.glm.logistic.hybrid" "intdata/56-GWAS_res_males_hardcall.PHENO1.glm.logistic.hybrid")

main_test=("ADD" "ADD")

mid_files=("intdata/56-gwas_females_hardcall_out.txt" "intdata/56-gwas_males_hardcall_out.txt")

sig_files=("intdata/56-gwas_females_hardcall_out_pfilt.txt" "intdata/56-gwas_males_hardcall_out_pfilt.txt")

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

## Filtering and Assembling Results

The previous section already did:

* Male-specific eur-only unrelated GWAS, 236 cases vs. 1180 controls. (5:1 matching)
* Female-specific eur-only unrelated GWAS, 473 cases vs. 946 controls. (2:1 matching)

And we applied MAF>0.01 filter to both cases and controls, limited to R2>0.8 SNPs, applied HWE p>1e-4 rule for controls, and rounded to 0/1/2 hardcalls.

Now we apply SFS-vs-CF filters (keep SNPs with p>1e-5 in the secondary GWAS), and apply pipeline filter (keep SNPs where chi-square test of AF difference between pipelines was >1e-4).

This is done in R:

```{r}
library(data.table)

## Control-Control results

con.file = paste0('intdata/47-gwas_out_noint_G_pfilt.txt')

conres = fread(con.file,data.table=FALSE)
colnames(conres)[1] = 'CHR'

cf_sfs_sug = conres[which(conres$P<1e-5),] # These are SNPs to remove

## Imputation Pipeline Results

imp.file = paste0('intdata/86-variants_with_differing_af_by_pipeline_CASES_OR_CONTROLS.txt')

impres = fread(imp.file,data.table=FALSE) 
colnames(impres)[1] = 'CHR'

imp_sug = impres[which(impres$chisq_p<1e-4),] # These are SNPs to remove

allrem = data.frame(v1=unique(c(cf_sfs_sug$ID,imp_sug$ID_batch)))
fwrite(allrem,file='91-snps_to_rem_concon_or_imppipeline.txt',col.names=FALSE)

########## Do everything for females

res.file = paste0('intdata/56-gwas_females_hardcall_out.txt')
sig.file = paste0('intdata/56-gwas_females_hardcall_out_pfilt.txt')

complete_allres = fread(res.file,data.table=FALSE)
complete_manres = fread(sig.file,data.table=FALSE)

final_allres = complete_allres[which(!(complete_allres$ID %in% cf_sfs_sug$ID | complete_allres$ID %in% imp_sug$ID_batch)),]

final_allres$ALT_FREQ = final_allres$A1_FREQ

final_allres$ALT_FREQ[which(final_allres$A1==final_allres$REF)]=1-final_allres$ALT_FREQ[which(final_allres$A1==final_allres$REF)]

fwrite(final_allres,'intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt',sep='\t')

########## Do everything for males

res.file = paste0('intdata/56-gwas_males_hardcall_out.txt')
sig.file = paste0('intdata/56-gwas_males_hardcall_out_pfilt.txt')

complete_allres = fread(res.file,data.table=FALSE)
complete_manres = fread(sig.file,data.table=FALSE)

final_allres = complete_allres[which(!(complete_allres$ID %in% cf_sfs_sug$ID | complete_allres$ID %in% imp_sug$ID_batch)),]

final_allres$ALT_FREQ = final_allres$A1_FREQ

final_allres$ALT_FREQ[which(final_allres$A1==final_allres$REF)]=1-final_allres$ALT_FREQ[which(final_allres$A1==final_allres$REF)]

fwrite(final_allres,'intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt',sep='\t')
```

## Mirrored Manhattan Plot

In R:

```{r}
library(hudson)
library(data.table)

male=fread('intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)
female=fread('intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)

topsnp = c('rs17318744','rs184652834','rs806795','rs4938014','rs4889844','rs2525035','rs197013')

# Use column names demanded by package:
colnames(male)[c(3,1,2,10)]=c("SNP", "CHR", "POS", "pvalue")
colnames(female)[c(3,1,2,10)]=c("SNP", "CHR", "POS", "pvalue")

male2 = male[which(male$pvalue <= 0.01 ),c("SNP", "CHR", "POS", "pvalue")]
female2 = female[which(female$pvalue <= 0.01 ),c("SNP", "CHR", "POS", "pvalue")]

my_chroms <- c(1:22, "X")

gmirror(top = male2, bottom = female2,
            tline=1e-5, bline=1e-5, 
            toptitle="JME GWAS in Males (n=1416)", bottomtitle = "JME GWAS in Females (n=1419)", 
            highlight_p = c(1e-5,1e-5), highlighter="green",
            chrcolor1 = "#AAAAAA",
            chrcolor2 = "#4D4D4D",
            chroms = my_chroms,
            file = "figures/my_miami_plot_refined")
```
