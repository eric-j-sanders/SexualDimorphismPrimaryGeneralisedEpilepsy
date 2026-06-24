# Purpose

When the codes in this folder are executed, the BIOJUME case data are merged with the SFS control data, and subset to those with inferred European ancestry. Then, relatedness is checked within the final merged dataset, to subset to unrelated individuals. Lastly, PCA is completed, the data are split to male- and female-specific subsets, and PCA-matching is completed to select the final GWAS sample populations.

# Input files

* The listed outputs of 01-CaseDataCleaning
* The listed outputs of 02-ControlDataCleaning
* Results from previous Grafpop ancestry inference on the BIOJUME sample saved to relative path intdata/10-BIOJUME_ancestry.csv

# Primary outputs

...

# Code

## Merge Cases and Controls

Using R, produce a list of IDs of European cases/controls:

```{r}
case = read.csv('intdata/10-BIOJUME_ancestry.csv')
control = read.csv('intdata/10-SFS2_ancestry.csv')

eur_ids = rbind(case[which(case$inferred_european==1),c(1,2)],control[which(control$inferred_european==1),c(1,2)])

write.table(eur_ids,file='intdata/10-eur_casecontrol_ids.txt',sep='\t',row.names=FALSE,quote=FALSE)
```
