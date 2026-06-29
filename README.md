# Project Overview

This repository contains copy of the code used to prepare data, complete analysis, and generate figures/tables as presented in the manuscript "Sexual Dimorphism in the Genetics of Primary Generalised Epilepsy" (link). This study completes sex-stratified GWAS of juvenile myoclonic epilepsy (JME), 

# Pipeline Description

The intention is not for the contents of this repo to allow out-of-the-box replication of results (see Data Availability statement below).

The original project had all code executed from the same directory with subfolders such as "rawdata" (for raw source data), "intdata" (for intermediate data), etc., but for clarity within this repo, the file structure has been altered to allow ease of review. This was done without changing the file pathing within the code itself, so the file paths in the code do not reflect the altered subfolders of this repo.

The codes are presented sequentially in subfolders:

1. Data cleaning: BIOJUME Cases
    * Case data filtered to individuals with probable/definite JME
    * formatted to facilitate merging
1. Data cleaning: SFS Controls
    * Control data has individuals removed with any psychiatric conditions
    * Ancestry inferred
    * formatted to facilitate merging
1. Data preparation (Merged case and control data)
    * Merge of case and control datasets
    * Subset to Europeans
    * Relatedness within entire merged set analyzed
    * PCA completed on unrelated European set
1. Secondary analyses to filter GWAS SNPs
    * Control-specific HWE filter
    * GWAS with alternative controls
    * GWAS with alternative imputation pipeline
1. Sex-Stratified GWAS
    * Sex-specific PC-based case-control matching
    * Sex-specific GWAS completed
    * Additional SNPs removed based on secondary analyses
    * Creating Manhattan plot
1. Male-Female Comparison of Summary Statistics
    * Cochran's Q Statistics and Heterogeneity Lambda GC
    * GPS-GEV Test
    * IVW Meta-Analysis
1. Gene-Based analysis
    * FUMA SNP2GENE gene annotation
    * FUMA GENE2FUNC enrichment analysis
1. Colocalization Analysis
    * Sex-Specific GWAS statistic subsets around peak loci
    * Sex-specific LD matrices around peak loci
    * Sex-specific eQTL statistics around peak loci
1. Polygenic Risk Score Calculation
    * Calculation of Neuroticism PRS in Case-Control cohort
    * Testing of association between neuroticism PRS and JME risk

# Data Access

(To copy statements from manuscript.)

# Computational Environment

The conda environment .yml file is included in this repository for referencing fine details of package versions.

Some primary programs and their versions:

* GNU bash, version 5.1.8(1)
* R version 3.6.3
* plink 1.90b7.7
* plink2 2.0.0a6.9 
* KING 2.2.7
* Grafpop version published 04/05/2021
* PCAmatchR version 0.3.3
* FUMA browser tool V2_r204_3510window