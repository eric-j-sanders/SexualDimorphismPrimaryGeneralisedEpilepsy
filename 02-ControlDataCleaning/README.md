# Purpose

When the codes in this folder are executed, the SFS control data are filtered to European-ancestry individuals without inferred/suspected/diagnosed psychiatric disorders, and have some formatting changes to facilitate merging with controls. Genetic ancestry is also inferred.

# Input files

* A set of pgen/psam/pvar files formatted for plink are assumed to be located with the relative path prefix "intdata/06-SFS2_prepared_genotype", representing the genotyped (not imputed) data from the SFS control population (before imputation). These data have already been through the QC process as detailed in the paper.

* A set of pgen/psam/pvar files formatted for plink are assumed to be located with the relative path prefix "intdata/06-SFS2_prepared_imputed", representing the post-imputation data from the SFS control population.

* A phenotype file with columns corresponding to survey entries of study participants, including reports of existing conditions or diagnoses, is assumed to be located with the relative path "rawdata/SFS2 - Data Pull - Epilepsy_Seizures_UPDATED December 23 2024.xlsx"

# Primary outputs

plink2 format files for both the genotyped and imputed data, that have been subset to individuals without neuropsychiatric conditions:

* pgen/psam/pvar files with relative pathing prefix intdata/07-SFS_prepared_genotype_nopsych
* pgen/psam/pvar files with relative pathing prefix intdata/07-SFS_prepared_imputed_nopsych
* A file listing inferred european ancestry status for each sample, with relative pathing intdata/10-SFS2_ancestry.csv

# Code

## Format Genotype+Imputed Data

In bash:

```{bash}
# Change psam files to have similar column structure for merging

awk 'BEGIN{FS=OFS="\t"} NR==1 {print "#FID", "IID", "SEX", "STATUS"} NR>1 {print "T" sprintf("%05d", NR-1), $0, 1}' intdata/06-SFS2_prepared_genotype.psam > SFS_with_FID.psam
mv SFS_with_FID.psam intdata/06-SFS2_prepared_genotype.psam

awk 'BEGIN{FS=OFS="\t"} NR==1 {print "#FID", "IID", "SEX", "STATUS"} NR>1 {print "T" sprintf("%05d", NR-1), $0, 1}' intdata/06-SFS2_prepared_imputed.psam > SFS_with_FID.psam
mv SFS_with_FID.psam intdata/06-SFS2_prepared_imputed.psam

# Remove header rows from pvar file

tail -n +45 intdata/06-SFS2_prepared_genotype.pvar | cut -f 1-5 > SFS_no_INFO.pvar
mv SFS_no_INFO.pvar intdata/06-SFS2_prepared_genotype.pvar

tail -n +43 intdata/06-SFS2_prepared_imputed.pvar | cut -f 1-5 > SFS_no_INFO.pvar
mv SFS_no_INFO.pvar intdata/06-SFS2_prepared_imputed.pvar
```

## Identify and Remove Reported Psychiatric Conditions

Use bash, convert xlsc doc to csv for easy use in R:

```{bash}
xlsx2csv --sheet 1 "rawdata/SFS2 - Data Pull - Epilepsy_Seizures_UPDATED December 23 2024.xlsx" | cut -d"," -f1,10-35 > intdata/06-SFS2_seiz_info.txt
```

And in R, read in survey responses, examine and assemble filters:

```{r}
library(dplyr)

pheno = read.csv('intdata/06-SFS2_seiz_info.txt',header=TRUE)
colnames(pheno) = c('IID','dx_other1','dx_other2','dx_seiz','seiz_fev1','seiz_fev0','seiz_fev_unknown','dx_epi','dx_BRE','dx_CAE','dx_JME','dx_WS','dx_LGS','dx_epi_other','dx_epi_unknown','seiz_focal_sp','seiz_focal_cp','seiz_focal_secgen','seiz_abs','seiz_gtc','seiz_gc','seiz_t','seiz_at','seiz_myo','seiz_unclass','seiz_unknown','dx_other3')
pheno$IID = toupper(pheno$IID)

fam = read.table('intdata/06-SFS2_prepared_genotype.psam',header=FALSE)

colnames(fam) = c("#FID",'IID','SEX','PHENO')

full = merge(fam,pheno,by='IID')

full = full[,c(2,1,3:ncol(full))]

# Produce some tables to see how the epilepsy info behaves/distributes
# These are manually observed to understand data distribution

combo_epi = full %>% count(dx_epi,dx_BRE,dx_CAE,dx_JME,dx_WS,dx_LGS,dx_epi_other,dx_epi_unknown) %>% arrange(desc(n))
print(combo_epi)

combo_seiz = full %>% count(dx_seiz,seiz_fev1,seiz_fev0,seiz_focal_sp,seiz_focal_cp,seiz_focal_secgen,seiz_abs,seiz_gtc,seiz_gc,seiz_t,seiz_at,seiz_myo) %>% arrange(desc(n))
print(combo_seiz)

combo_major = full %>% count(dx_seiz,seiz_fev1,seiz_fev0,dx_epi,dx_epi_unknown) %>% arrange(desc(n))
print(combo_major)

# It is deemed sufficient to just remove anyone who reports seizures/epilepsy in the main two variables

full = full[which(is.na(full$dx_seiz) | full$dx_seiz==0),]

# Also check for neuropsychiatric conditions:
others=full[which(full$dx_other1!="" | full$dx_other2!="" | full$dx_other3!=""),c(2,5,6,30)]

# Manually assemble list of character strings that are used to mark individuals with neuropsychiatric conditions
exclusions = c('sensory','autism','adhd','schizo','anxiety','ocd','memory','cognitive','frontal lobe','spd','asd','hallucinations','bipolar','down syndrome','psychosis','delay','mosaic turner','odd','brain injury','chromosoem 15','appxia','apraxia','dyspraxia', 'mutism','dsylexia','dislexia','dyslexia','audio processing','auditory processing','visual perception','did','cerebral palsy', 'ptsd','adjustment disorder', 'aphantasia','developmental coordination disorder','developmental coordination  disorder','cognitive affective disorder','hydrocephalus','Atenshen disorder','tdah','tourettes','auditory processing')

# Two samples report the "asd" string where context suggests it does not refer to autism spectrum disorder --> Ensure those are not removed
full$dx_other3[which(full$dx_other3=='ASD  VSD  HRHS')] = "Atrial septal defect VSD HRHS" # DONT remove ASD VSD HRHS patient, ASD does not refer to autism in this case
full$dx_other3[which(full$dx_other3=='asd heart ')] = "Atrial septal defect heart" # DONT remove ASD VSD HRHS patient, ASD does not refer to autism in this case

# Remove all psychiatric conditions now:
full_reduced = full[sapply(full$dx_other1, function(x) !any(grepl(paste(exclusions,collapse='|'), x, ignore.case=TRUE))),]
full_reduced = full_reduced[sapply(full_reduced$dx_other3, function(x) !any(grepl(paste(exclusions,collapse='|'), x, ignore.case=TRUE))),]

# Make a file to quickly glance at the diagnoses that have NOT been filtered at this point:
others_reduced=full_reduced[which(full_reduced$dx_other1!="" | full_reduced$dx_other2!="" | full_reduced$dx_other3!=""),c(2,5,6,30)]
# Manual observation confirms no remaining survey results suggest anyone with neuropsychiatric conditions.

kept = unique(c(others_reduced$dx_other1,others_reduced$dx_other2,others_reduced$dx_other3))

write.table(kept,file='intdata/06-SFS2_diagnoses_retained.txt',quote=FALSE,row.names=FALSE)

write.table(full_reduced[,1:4],'intdata/06-SFS2_psych_removed.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")
```

Now in bash these samples are removed from both the genotype and imputed data sets:

```{bash}
plink2 --pfile intdata/06-SFS2_prepared_genotype --keep intdata/06-SFS2_psych_removed.txt --make-pgen --out intdata/07-SFS2_prepared_genotype_nopsych

plink2 --pfile intdata/06-SFS2_prepared_imputed --keep intdata/06-SFS2_psych_removed.txt --make-pgen --out intdata/07-SFS2_prepared_imputed_nopsych
```

# Ancestry analysis

```{bash}
plink2 --pfile intdata/07-SFS2_prepared_genotype_nopsych --make-bed --out intdata/07-forgrafpop

grafpop intdata/07-forgrafpop.bed intdata/10-SFS2-grafpop.txt

perl ~/Programs/PlotGrafPopResults.pl intdata/10-SFS2-grafpop.txt figures/SFS2_GRAFPOP_output.png

perl ~/Programs/SaveSamples.pl intdata/10-SFS2-grafpop.txt intdata/10-SFS2-grafpop_samples.txt

rm intdata/07-forgrafpop*
```

And, using R, assembling a tidy dataset of ancestry analysis results for SFS:

```{r}
psam = read.table('intdata/07-SFS2_prepared_genotype_nopsych.psam')

grafpop = read.table('intdata/10-SFS2-grafpop_samples.txt',comment.char='',header=TRUE,sep='\t')

merged = merge(psam,grafpop,by.x='V2',by.y='Subject')

merged$inferred_european = as.numeric(merged$Computed.population=='European')

merged = merged[,c('V1','V2','PopID','Computed.population','inferred_european')]

colnames(merged)[1:2]=c('FID',"IID")

write.csv(merged,file='intdata/10-SFS2_ancestry.csv',col.names=TRUE,row.names=FALSE)
```