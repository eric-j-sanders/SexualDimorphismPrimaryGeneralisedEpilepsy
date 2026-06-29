# Purpose

When the codes in this folder are executed, data are prepared to upload to FUMA's SNP2GENE and GENE2FUNC pipeline for enrichment analysis.

# Input files

* The listed outputs of 05-SexStratifiedGWAS

# Primary outputs

* Gene annotation results saved to FUMA_Output/V2_r204_3510window/{male/female}_GWAS_SNP2GENE/genes.txt
* Enrichment analysis results saved to FUMA_Output/V2_r204_3510window/{male/female}_GWAS_gene2func/GS.txt

# Code

## Cleaning for FUMA

We need a sample size column and for the summary statistics to be zipped:

```{bash}
awk -F'\t' 'BEGIN {OFS="\t"} {print (NR==1 ? "SNP" : $3), $6, $10, (NR==1 ? "N" : "1419")}' intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt > intdata/95-females_forMAGMA.txt

bgzip intdata/95-females_forMAGMA.txt

awk -F'\t' 'BEGIN {OFS="\t"} {print (NR==1 ? "SNP" : $3), $6, $10, (NR==1 ? "N" : "1416")}' intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt > intdata/95-males_forMAGMA.txt

bgzip intdata/95-males_forMAGMA.txt
```

## Upload to SNP2GENE

These files are uploaded to FUMA SNP2GENE with the following settings: 

* In Section 1 the relevant gzip file is uploaded (created in the previous code block) 
* In section 2 the Sample Size is specified via column name “N”, and the “Maximum P-value of lead SNPs” is increased up to 1e-5, and the r2 threshold to define independent significant SNPs is reduced to 0.4 
* In section 3-2 select to perform eQTL mapping, selecting tissue types for GTEx v8 Brain and Nerve tissues, as well as GTEx v8 Whole Blood (Whole Blood in GTEx known to be more powerful from higher sample size) 
* In section 3-3 select chromatin interaction mapping, and select the chromatin interaction data for HiC brain tissues (adult cortex, fetal cortex, prefrontal cortex, hippocampus, neural progenitor cell). Under “annotate enhance/promoter regions” select all Brain tissues, then check the boxes undereneath “Filter SNPs by enhancers” and “filter SNPs by promoters” 
* In Section 4 the Ensembl version is changed to v110 and the Gene type category selects Protein Coding and Processed transcripts and lncrna (not pseudogenes, IG genes, etc.) 
* In section 6 “Perform MAGMA” is selected, with the Gene windows set to “35,10”, and with both GTEx v8 options and with Brainspan set of 29 age options

Results are saved to FUMA_Output/V2_r204_3510window/{male/female}_GWAS_SNP2GENE/genes.txt

## Upload to GENE2FUNC

The results of the previous section are piped into FUMA GENE2FUNC, with these setting changes from default: 

* Background genes limited to Protein coding, lncRNA, and Processed Transcripts. 
* Ensemble version 110 
* Gene expression data sets include GTEx v8 (both options) and Brainspan 29 different ages. 
* MHC region exclusion box is checked.

Results are saved to FUMA_Output/V2_r204_3510window/{male/female}_GWAS_gene2func/GS.txt