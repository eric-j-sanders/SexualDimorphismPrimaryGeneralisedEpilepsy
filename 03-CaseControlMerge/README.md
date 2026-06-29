# Purpose

When the codes in this folder are executed, the BIOJUME case data are merged with the SFS control data, and subset to those with inferred European ancestry. Then, relatedness is checked within the final merged dataset, to subset to unrelated individuals. Lastly, PCA is completed, the data are split to male- and female-specific subsets, and PCA-matching is completed to select the final GWAS sample populations.

# Input files

* The listed outputs of 01-CaseDataCleaning
* The listed outputs of 02-ControlDataCleaning
* Results from previous Grafpop ancestry inference on the BIOJUME sample saved to relative path intdata/10-BIOJUME_ancestry.csv
* The shell scripts saved to this directory
* An "extended_ld_regions.txt" file that lists extended LD regions (HG38) for removal before completing PCA.

# Primary outputs

* A set of pgen/psam/pvar files in plink format representing all the European-ancestry merged case and control imputed genotype data, for the JME-positive patients and controls without neuropsychiatric conditions, with relative pathing prefix intdata/14-merge_imputed_eur_unrel
* A set of pgen/psam/pvar files in plink format representing all the European-ancestry merged case and control genotyped data (not imputed), for the JME-positive patients and controls without neuropsychiatric conditions, with relative pathing prefix intdata/14-merge_genotype_eur_unrel
* PCA results on the noted merged population, saved to relative paths intdata/16-pca_res.eigenvec and intdata/16-pca_res.eigenval
* A text file containing the sex, case/control status, and first 10 PCs for all the merged cases and controls, saved to relative path intdata/17-eur_unrel_pheno.txt

# Code

## Merge Cases and Controls

Using R, produce a list of IDs of European cases/controls:

```{r}
case = read.csv('intdata/10-BIOJUME_ancestry.csv')
control = read.csv('intdata/10-SFS2_ancestry.csv')

eur_ids = rbind(case[which(case$inferred_european==1),c(1,2)],control[which(control$inferred_european==1),c(1,2)])

write.table(eur_ids,file='intdata/10-eur_casecontrol_ids.txt',sep='\t',row.names=FALSE,quote=FALSE)
```

### Merge Genotype Data

Merge cases and controls:

```{bash}
sbatch code/merge_bplink_files_v2_for_plink2.sh intdata/07-SFS2_prepared_genotype_nopsych intdata/09-BIOJUME_prepared_genotype intdata/11-merge_genotype genmerge
```

Use R to tidy the genotype psam post-merge:

```{r}
psam = read.table('intdata/11-merge_genotype.psam',comment.char='',header=TRUE)
psam$`#FID` = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][1])
psam$IID = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][2])

psam$`#FID` = gsub('-','_',psam$`#FID`)
psam$IID = gsub('-','_',psam$IID)

psam2 = psam[,c('#FID','IID','SEX')]

psam_bj = read.table('intdata/09-BIOJUME_prepared_genotype.psam',comment.char='',header=TRUE)
colnames(psam_bj)[1]='#FID'

psam_sfs = read.table('intdata/07-SFS2_prepared_genotype_nopsych.psam',comment.char='',header=TRUE)
colnames(psam_sfs)[1]='#FID'

psam_ogs = rbind(psam_bj,psam_sfs)

psam2$STATUS=NA

psam2$SEX[match(psam_ogs$IID, psam2$IID)] = psam_ogs$SEX
psam2$STATUS[match(psam_ogs$IID, psam2$IID)] = psam_ogs$STATUS

write.table(psam2,'intdata/11-merge_genotype.psam',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
```

And subset to European samples using bash:

```{bash}
plink2 --pfile intdata/11-merge_genotype --keep intdata/10-eur_casecontrol_ids.txt --make-pgen --out intdata/12-merge_genotype_eur
```

### Merge Imputed Data

In bash:

```{bash}
sbatch code/merge_bplink_files_v2_for_plink2.sh intdata/07-SFS2_prepared_imputed_nopsych intdata/09-BIOJUME_prepared_imputed intdata/11-merge_imputed impmerge
```

And tidying in R:

```{r}
psam = read.table('intdata/11-merge_imputed.psam',comment.char='',header=TRUE)
psam$`#FID` = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][1])
psam$IID = sapply(psam$X.IID,function(x) strsplit(x,'_')[[1]][2])

psam$`#FID` = gsub('-','_',psam$`#FID`)
psam$IID = gsub('-','_',psam$IID)

psam2 = psam[,c('#FID','IID','SEX')]

psam_bj = read.table('intdata/09-BIOJUME_prepared_imputed.psam',comment.char='',header=TRUE)
colnames(psam_bj)[1]='#FID'

psam_sfs = read.table('intdata/07-SFS2_prepared_imputed_nopsych.psam',comment.char='',header=TRUE)
colnames(psam_sfs)[1]='#FID'

psam_ogs = rbind(psam_bj,psam_sfs)

psam2$STATUS=NA

psam2$SEX[match(psam_ogs$IID, psam2$IID)] = psam_ogs$SEX
psam2$STATUS[match(psam_ogs$IID, psam2$IID)] = psam_ogs$STATUS

write.table(psam2,'intdata/11-merge_imputed.psam',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
```

And subsetting to Europeans:

```{bash}
plink2 --pfile intdata/11-merge_imputed --keep intdata/10-eur_casecontrol_ids.txt --make-pgen --out intdata/12-merge_imputed_eur
```

## Relatedness in Merged Population

Relatedness is checked using king in bash:

```{sh}
plink2 --pfile intdata/12-merge_genotype_eur --make-bed --out intdata/13-forking
plink --bfile intdata/13-forking --make-bed --out intdata/13-forking2

king -b intdata/13-forking2.bed --related --degree 2 --prefix intdata/13-kinship_all
king -b intdata/13-forking2.bed --unrelated --degree 2 --prefix intdata/13-kinship_unrel

rm intdata/13-forking*

plink2 --pfile intdata/12-merge_genotype_eur --keep intdata/13-kinship_unrelunrelated.txt --make-pgen --out intdata/14-merge_genotype_eur_unrel
```

Imputed data is also subset:

```{bash}
plink2 --pfile intdata/12-merge_imputed_eur --keep intdata/13-kinship_unrelunrelated.txt --make-pgen --out intdata/14-merge_imputed_eur_unrel
```

## PCA in Merged Population

In bash:

```{bash}
bash code/prePCA_plink2.sh intdata/14-merge_genotype_eur_unrel intdata/15

bash code/plink2_PCA.sh 'intdata/15-pruned' 'intdata/16-pca_res'
```

PCA was done on 58k SNPs.

Use R to save a tidy table of first 10 PCs with sex and phenotype columns:

```{r}
psam = read.table('intdata/14-merge_genotype_eur_unrel.psam')
colnames(psam)=c('FID','IID','SEX','STATUS')

eigenval = read.table('intdata/16-pca_res.eigenval')

eigenvec = read.table('intdata/16-pca_res.eigenvec')
colnames(eigenvec) = c('FID','IID',paste0('PC',1:(ncol(eigenvec)-2)))

outfile = 'intdata/17-eur_unrel_pheno.txt'

eigenval$V2 = 1:nrow(eigenval)
eigenval$V3 = eigenval$V1/sum(eigenval$V1)

new_pheno = merge(psam, eigenvec, by = "IID", all = FALSE)[,c(2,1,3,4,6:15)]
colnames(new_pheno)[1] = 'FID'

new_pheno$Pheno = factor(new_pheno$STATUS,levels=c(1,2),labels=c('SFS','BIOJUME'))

new_pheno$IID = gsub('-','_',new_pheno$IID)
new_pheno$FID = gsub('-','_',new_pheno$FID)

write.table(new_pheno[,1:14], outfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
```