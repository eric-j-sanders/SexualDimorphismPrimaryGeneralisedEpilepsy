# Purpose

When the codes in this folder are executed, neuroticism polygenic risk scores are calculated in the case-control populations, and the PRS measures are compared to JME risk via regression analysis.

# Input files

* The listed outputs of 05-SexStratifiedGWAS
* The listed outputs of 04-SecondaryFilterCollection
* The neuroticism GWAS statistics for PRS calculation, available in rawdata/11-gwas_stats_wor_final.txt (worry) and rawdata/11-gwas_stats_dep_final.txt (depressed affect).

# Primary outputs

* Logistic regression models testing association between depressed affect/worry neuroticism and JME risk in males/females (four models total).

# Code

## Data Preparation

First take the sex-stratified pre-GWAS data and remove the SNPs marked from the Secondary Filter Collection step:

In R:

```{r}
library(data.table)

cf_sfs_sug = fread('intdata/47-gwas_out_noint_G_pfilt.txt',data.table=FALSE)
cf_sfs_sug = cf_sfs_sug[which(cf_sfs_sug$P<1e-5),]

impsug = fread('intdata/86-variants_with_differing_af_by_pipeline_CASES_OR_CONTROLS.txt',data.table=FALSE)
impsug = impsug[which(impsug$chisq_p<1e-4),]

removal_ids = c(cf_sfs_sug$ID,impsug$ID_batch)

fwrite(data.frame(id=removal_ids),'intdata/97-cfsfs_and_imp_snps_toremove.txt',sep='\t',col.names=FALSE)
```

Then in bash:

```{bash}
plink2 --bfile intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls --exclude intdata/97-cfsfs_and_imp_snps_toremove.txt --make-bed --out intdata/97-males_allfilters

plink2 --bfile intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls --exclude intdata/97-cfsfs_and_imp_snps_toremove.txt --make-bed --out intdata/97-females_allfilters
```

Filter to SNPs in common between neuroticism GWAS and case-control sample:

First in R:

```{r}
library(data.table)
library(dplyr)

for(type in c('dep','wor')){
    GWAS = fread(paste0("rawdata/11-gwas_stats_",type,"_final.txt"), header=TRUE,data.table=FALSE)

    names(GWAS)[2] = "BP_37"
    names(GWAS)[15] = "BP"
    for(sex in c('male','female')){
        bim = fread(paste0("intdata/97-",sex,"s_allfilters.bim"), header=FALSE,data.table=FALSE)
        snps = bim[,"V2"]
        
        GWAS_Final = GWAS[which(GWAS$RSID %in% snps),]
        
        bim = bim[which(bim$V2 %in% GWAS_Final$RSID),]
        
        fwrite(GWAS_Final, paste0("intdata/98-neurgwas_stats_",type,"_final_common_",sex,".txt"), sep="\t", row.names = FALSE, quote = FALSE)
        
        SNP_P = GWAS_Final[,c(3,12)]

        fwrite(SNP_P, paste0("intdata/98-neurgwas_stats_",type,"_final_common_",sex,"_P.txt"), sep="\t", row.names = FALSE, quote = FALSE)
        
        snp_extraction_file = data.frame(CHR=GWAS_Final$CHR, START=GWAS_Final$BP, END=GWAS_Final$BP,NAME=GWAS_Final$RSID)
  
    fwrite(snp_extraction_file, paste0('intdata/98-snp_extraction_',type,'_',sex,'.txt'),quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
    }
}
```

Then in bash:

```{bash}
diff -qs intdata/98-snp_extraction_dep_male.txt intdata/98-snp_extraction_wor_male.txt # Files don't change by dep/wor
diff -qs intdata/98-snp_extraction_dep_female.txt intdata/98-snp_extraction_wor_female.txt # Files don't change by dep/wor
diff -qs intdata/98-snp_extraction_dep_male.txt intdata/98-snp_extraction_dep_female.txt # Files change by sex

# Therefore only need to do SNP extraction with "dep" files, the "wor" extraction would be identical

plink2 --bfile intdata/97-females_allfilters --extract range intdata/98-snp_extraction_dep_female.txt --make-bed --out intdata/98-females_matchsnps

plink2 --bfile intdata/97-males_allfilters --extract range intdata/98-snp_extraction_dep_male.txt --make-bed --out intdata/98-males_matchsnps


```

## PRS Calculation

As in Liuhanen et al, the PRS is calculated using PRSice program, using the LD-clumping commands (−-clump-kb 500, −-clump-p 1.000000 and –clump-r2 0.250000) and the p-value threshold of 0.1.

```{bash}
# Male dep
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base intdata/98-neurgwas_stats_dep_final_common_male.txt \
    --target intdata/98-males_matchsnps \
    --clump-kb 500 \
    --clump-p 1.000000 \
    --clump-r2 0.250000 \
    --pvalue P \
    --bar-levels 0.1 --fastscore --no-full --all-score --no-regress \
    --stat BETA \
    --snp RSID \
    --chr CHR \
    --bp BP \
    --a1 A1 \
    --a2 A2 \
    --out intdata/99-PRS_output_male_dep
    
# Male wor
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base intdata/98-neurgwas_stats_wor_final_common_male.txt \
    --target intdata/98-males_matchsnps \
    --clump-kb 500 \
    --clump-p 1.000000 \
    --clump-r2 0.250000 \
    --pvalue P \
    --bar-levels 0.1 --fastscore --no-full --all-score --no-regress \
    --stat BETA \
    --snp RSID \
    --chr CHR \
    --bp BP \
    --a1 A1 \
    --a2 A2 \
    --out intdata/99-PRS_output_male_wor   
    
# Female dep
Rscript PRSice.R \ 
    --prsice PRSice_linux \
    --base intdata/98-neurgwas_stats_dep_final_common_female.txt \
    --target intdata/98-females_matchsnps \
    --clump-kb 500 \
    --clump-p 1.000000 \
    --clump-r2 0.250000 \
    --pvalue P \
    --bar-levels 0.1 --fastscore --no-full --all-score --no-regress \
    --stat BETA \
    --snp RSID \
    --chr CHR \
    --bp BP \
    --a1 A1 \
    --a2 A2 \
    --out intdata/99-PRS_output_female_dep
    
# Female wor
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base intdata/98-neurgwas_stats_wor_final_common_female.txt \
    --target intdata/98-females_matchsnps \
    --clump-kb 500 \
    --clump-p 1.000000 \
    --clump-r2 0.250000 \
    --pvalue P \
    --bar-levels 0.1 --fastscore --no-full --all-score --no-regress \
    --stat BETA \
    --snp RSID \
    --chr CHR \
    --bp BP \
    --a1 A1 \
    --a2 A2 \
    --out intdata/99-PRS_output_female_wor   
```

## Testing for PRS association with JME risk

In R:

```{r}
library(data.table)
library(ggplot2)
library(dplyr)

md = fread('intdata/99-PRS_output_male_dep.all_score')
colnames(md)[3]='PRS_dep'
mw = fread('intdata/99-PRS_output_male_wor.all_score')
colnames(mw)[3]='PRS_wor'
fd = fread('intdata/99-PRS_output_female_dep.all_score')
colnames(fd)[3]='PRS_dep'
fw = fread('intdata/99-PRS_output_female_wor.all_score')
colnames(fw)[3]='PRS_wor'

prs_m = merge(md,mw,by=c('FID','IID'))
prs_m$FID = gsub('_','-',prs_m$FID)
prs_m$IID = gsub('_','-',prs_m$IID)
prs_f = merge(fd,fw,by=c('FID','IID'))
prs_f$FID = gsub('_','-',prs_f$FID)
prs_f$IID = gsub('_','-',prs_f$IID)

mfam = fread('intdata/97-males_allfilters.fam')[,c(1,2,5,6)]
colnames(mfam)=c('FID','IID','Sex','JME')
mfam$FID = gsub('_','-',mfam$FID)
mfam$IID = gsub('_','-',mfam$IID)
ffam = fread('intdata/97-females_allfilters.fam')[,c(1,2,5,6)]
colnames(ffam)=c('FID','IID','Sex','JME')
ffam$FID = gsub('_','-',ffam$FID)
ffam$IID = gsub('_','-',ffam$IID)

mcov = fread('intdata/51-cov_matched_males_10.txt')
mcov$FID = gsub('_','-',mcov$FID)
mcov$IID = gsub('_','-',mcov$IID)
fcov = fread('intdata/51-cov_matched_females_10.txt')
fcov$FID = gsub('_','-',fcov$FID)
fcov$IID = gsub('_','-',fcov$IID)

mcov_merge = merge(mcov,mfam,by=c('FID','IID'))
fcov_merge = merge(fcov,ffam,by=c('FID','IID'))

mfin = merge(mcov_merge,prs_m,by=c('FID','IID'))
ffin = merge(fcov_merge,prs_f,by=c('FID','IID'))

mfin$PRS_dep = (mfin$PRS_dep-mean(mfin$PRS_dep))/sd(mfin$PRS_dep)
mfin$PRS_wor = (mfin$PRS_wor-mean(mfin$PRS_wor))/sd(mfin$PRS_wor)
mfin$Status = factor(mfin$JME,levels=c(1,2),labels=c('Control (Spit2 Group)','Case (BIOJUME JME)'))
mfin$Sex = factor(mfin$Sex, levels=c(1,2),labels=c('Male','Female'))

ffin$PRS_dep = (ffin$PRS_dep-mean(ffin$PRS_dep))/sd(ffin$PRS_dep)
ffin$PRS_wor = (ffin$PRS_wor-mean(ffin$PRS_wor))/sd(ffin$PRS_wor)
ffin$Status = factor(ffin$JME,levels=c(1,2),labels=c('Control (Spit2 Group)','Case (BIOJUME JME)'))
ffin$Sex = factor(ffin$Sex, levels=c(1,2),labels=c('Male','Female'))

summary(glm(Status~PRS_dep+PC1+PC2+PC3+PC4+PC5,data=mfin,family=binomial(link='logit')))
summary(glm(Status~PRS_wor+PC1+PC2+PC3+PC4+PC5,data=mfin,family=binomial(link='logit')))
summary(glm(Status~PRS_dep+PC1+PC2+PC3+PC4+PC5,data=ffin,family=binomial(link='logit')))
summary(glm(Status~PRS_wor+PC1+PC2+PC3+PC4+PC5,data=ffin,family=binomial(link='logit')))
```