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
    * PCs used to match cases and controls
1. Secondary analyses to filter GWAS SNPs
    * Control-specific HWE filter
    * GWAS with alternative controls
    * GWAS with alternative imputation pipeline
1. Sex-Stratified GWAS




# Data Access

# Computational Environment