# Purpose

When the codes in this folder are executed, three secondary analytical steps are completed to prepare additional SNP filters for application to the GWAS results:

* HWE testing specifically in controls (SNPs will be removed if they have P<1e-4 in this test)
* A control vs. control GWAS where both groups are of European ancestry (SNPs will be removed if they have P<1e-5 in this GWAS)
* Chi-square testing of AF differences between completing imputation separately by batch or in merged batches (SNPs will be removed if both imputation pipelines have imputation R2>0.8 but the SNP has P<1e-4 in the chi-square test of AF difference)

# Input files

* The listed outputs of all previous repository blocks
* CF input...
* Imputation input...
* The shell scripts saved to this directory
* An "extended_ld_regions.txt" file that lists extended LD regions (HG38) for removal before completing PCA.

# Primary outputs

* 

# Code

## HWE Testing in Controls

Copy from Section 17

## Secondary GWAS Between SFS and CF Controls

Copy from Section 23

## Comparison of Imputation Pipeline (Identifying SNPs with significant AF variation)

Copy from Sections 34, 37