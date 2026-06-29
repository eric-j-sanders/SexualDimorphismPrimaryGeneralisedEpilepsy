# Purpose

When the codes in this folder are executed, data are organized to complete sex-specific colocalization testing for all 9 GWAS peaks of interest, using sex-specific GTEx v10 data and in-sample sex-specific LD matrices.

# Input files

* The listed outputs of 05-SexStratifiedGWAS
* Sex-Specific eQTL data, already subset to one file per testing region + gene + tissue combination. For one testing sex, region, gene, tissue, there will be a list of SNPs with corresponding eQTL test statistics. The male files are saved to intdata/102-GTEx_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/\*.txt.gz. The female files are saved to intdata/102-GTEx_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/\*.txt.gz
* Code scripts listed in this folder

# Primary outputs

For each sex, for each of the 9 GWAS peaks highlighted in the manuscript, a trio of files are created to upload to the LocusFocus web tool for colocalization:

* Sex-specific eQTL statistics: intdata/94-${sex}_${reg}_snps_eqtlv10.html
* Sex-specific LD matrix: intdata/94-${sex}_${reg}_snps_ldcorrected.ld
* Sex-specific GWAS statistics: intdata/94-${sex}_${reg}_res.txt

# Code

## Subsetting GWAS Statistics to Regions of Interest:

In R:

```{bash}
regions=data.frame(
    rsid=c('rs17318744','rs184652834','rs806795','rs4938014','rs4889844','rs2525035','rs197013','rs702041','rs8058223'),
    chr=c(1,3,6,11,17,22,'X',3,16),
    pos=c(210071125,14941254,26205065,113393799,80325697,19094833,28650816,124717498,15810351),
    sex=c('m','f','m','f','f','m','f','f','m')
)

regions$start = regions$pos-200000
regions$end = regions$pos+200000

write.table(regions,file='intdata/92-regions_of_interest.txt', sep = "\t", row.names = FALSE, quote = FALSE, eol = "\n",col.names=FALSE) 

write.table(regions[,c(2,5,6)],file='intdata/92-regions_of_interest_for_plink_selection.txt', sep = "\t", row.names = FALSE, quote = FALSE, eol = "\n",col.names=FALSE) 
```

And in bash:

```{bash}
### FEMALES

in_regions="intdata/92-regions_of_interest_for_plink_selection.txt"
in_gwas="intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt"

while read -r chr start end; do
  out1="intdata/93-female_${chr}_${start}_${end}_snps.txt"
  out2="intdata/94-female_${chr}_${start}_${end}_res.txt"
  
  awk -v chr="$chr" -v start="$start" -v end="$end" '
    NR==1 { next }
    $1==chr && $2>=start && $2<=end { print $3 }
  ' "$in_gwas" > "$out1"
  
  awk -v chr="$chr" -v start="$start" -v end="$end" '
    NR==1 { print $0; next }
    $1==chr && $2>=start && $2<=end { print $0 }
  ' "$in_gwas" > "$out2"

done < "$in_regions"

for f in intdata/93-female*_snps.txt; do
    base="${f%.txt}"
    plink2 \
        --bfile intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
        --extract "$f" \
        --make-bed \
        --out "$base"
    echo "........................................"
done

##### Female peaks, female sample:
# Chr3 region1: 1062 SNPs
# Chr3 region2: 1557 SNPs
# Chr11 region: 908 SNPs
# Chr17 region: 1139 SNPs
# ChrX region: 519 SNPs

##### Male peaks, female sample:
# Chr1 region: 630 SNPs
# Chr6 region: 1248 SNPs
# Chr16 region: 1002 SNPs
# Chr22 region: 605 SNPs


### MALES

in_regions="intdata/92-regions_of_interest_for_plink_selection.txt"
in_gwas="intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt"

while read -r chr start end; do
  out1="intdata/93-male_${chr}_${start}_${end}_snps.txt"
  out2="intdata/94-male_${chr}_${start}_${end}_res.txt"
  
  awk -v chr="$chr" -v start="$start" -v end="$end" '
    NR==1 { next }
    $1==chr && $2>=start && $2<=end { print $3 }
  ' "$in_gwas" > "$out1"

    awk -v chr="$chr" -v start="$start" -v end="$end" '
    NR==1 { print $0; next }
    $1==chr && $2>=start && $2<=end { print $0 }
  ' "$in_gwas" > "$out2"

done < "$in_regions"

for f in intdata/93-male*_snps.txt; do
    base="${f%.txt}"
    plink2 \
        --bfile intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
        --extract "$f" \
        --make-bed \
        --out "$base"
    echo "........................................"
done

##### Female peaks, female sample:
# Chr3 region: 1059 SNPs
# Chr3 region2: 1503 SNPs
# Chr11 region: 898 SNPs
# Chr17 region: 1139 SNPs
# ChrX region: 484 SNPs

##### Male peaks, male sample
# Chr1 region: 668 SNPs
# Chr6 region: 1265 SNPs
# Chr16 region: 1021 SNPs
# Chr22 region: 603 SNPs
```

## Calculate Sex-Specific LD Matrices, Apply Corrections

In bash:

```{bash}
for prefix in $(ls intdata/93-*_snps.bed | sed 's/\..*//' | sort -u); do
    plink2 --bfile "$prefix" --r-unphased square --out "$prefix" 
done
```

These LD matrices are in need of regularization/correction to address small sample bias and any numerical issues. The steps will be:

* Apply Bulik-Sullivan correction replace r2 with r2-(1-r2)/(n-2) as in Bulik-Sullivan, B. K. et al. LD Score regression distinguishes confounding from polygenicity in genomewide association studies. Nature Genetics 47, 291–295 (2015).
* Ensure diagonal is fixed at 1 and symmetry is forced
* Employ eigenvalue regularization by enforcing that the smallest eigenvalues are no lower than 1e-4

This can be done in R:

```{r}
library(data.table)
options(scipen=999) #Don't want scientific notation in saved LD matrix

all_files=list.files("intdata/",pattern="93-")
prefixes = unique(sapply(strsplit(all_files,'_snps'),'[[',1))

for(pre in prefixes){
    # Load in LD matrix
    print(paste0('Starting ',pre))
    
    ld = fread(paste0('intdata/',pre,'_snps.unphased.vcor1'),data.table=FALSE)
    
    ld = as.matrix(ld)
    
    # Employ Bulik-Sullivan correction
    r2 = ld^2
    if(grepl('female',pre)){
        n=1419
    }else{
        n=1416
    }
    
    r2_corr = r2-(1-r2)/(n-2)
    
    r2_corr = pmax(r2_corr,0)
    
    # Force proper symmetry/diagonal
    
    diag(r2_corr) = 1
    r2_corr = (r2_corr+t(r2_corr))/2
    
    # Eigenvalue regularization
    
    ed = eigen(r2_corr, symmetric=TRUE)
    ed$values = pmax(ed$values, 1e-4)
    r2_reg = ed$vectors %*% diag(ed$values) %*% t(ed$vectors)
    
    d_scale = 1/sqrt(diag(r2_reg))
    
    r2_reg = diag(d_scale) %*% r2_reg %*% diag(d_scale)
    
    r2_reg = pmin(pmax(r2_reg, 0), 1)
    diag(r2_reg)=1
    r2_reg = (r2_reg+t(r2_reg))/2
    
    r2_reg = round(r2_reg,15)
    
    newpre = sub('93-','94-',pre)
    out=paste0('intdata/',newpre,'_snps_ldcorrected.ld')
    print(paste0("... Saving ",out))
    fwrite(r2_reg,file=out,sep='\t',eol='\n',col.names=FALSE)
}
```

## Processing Sex-Specific eQTL data in Regions of interest

Prepare "descriptor files" that list the datasets corresponding to a particular region (a list of files, where each file represents one gene-tissue pair tested in the region):

```{bash}
# Loop through all the unique regions and collect the corresponding files based on chr_start_end, making a different descriptor file per region.

while IFS=$'\t' read -r col1 col2 col3; do
    namestring="${col1}_${col2}_${col3}"
    echo "$namestring" 
    
    ls intdata/102-GTExv10_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/*"${namestring}"*.txt.gz > intdata/102-GTExv10_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt
    
    # Finish descriptor files:

        awk -F'\t' 'BEGIN {OFS="\t"} {
          n = split($1, path_parts, "/");
          base = path_parts[n];
          # Remove everything from the 3rd underscore-separated chunk onward (coords + ext)
          # Instead, split on "_ENSG" to isolate tissue and gene separately
          split(base, halves, "_ENSG");
          tissue = halves[1];
          # Extract gene id: take ENSG part, then grab only "ENSGxxxxxxx.xx" (stop before next underscore)
          split("ENSG" halves[2], gene_parts, "_");
          gene = gene_parts[1];
          col2 = tissue "_" gene;
          print $1, col2, "chr", "pos", "variant_id", "pval_nominal"
        }' intdata/102-GTExv10_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt > temp_file.txt && mv temp_file.txt intdata/102-GTExv10_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt

done < "intdata/92-regions_of_interest_for_plink_selection.txt"

while IFS=$'\t' read -r col1 col2 col3; do
    namestring="${col1}_${col2}_${col3}"
    echo "$namestring" 
    
    ls intdata/102-GTExv10_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/*"${namestring}"*.txt.gz > intdata/102-GTExv10_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt
    
    # Finish descriptor files:

        awk -F'\t' 'BEGIN {OFS="\t"} {
          n = split($1, path_parts, "/");
          base = path_parts[n];
          # Remove everything from the 3rd underscore-separated chunk onward (coords + ext)
          # Instead, split on "_ENSG" to isolate tissue and gene separately
          split(base, halves, "_ENSG");
          tissue = halves[1];
          # Extract gene id: take ENSG part, then grab only "ENSGxxxxxxx.xx" (stop before next underscore)
          split("ENSG" halves[2], gene_parts, "_");
          gene = gene_parts[1];
          col2 = tissue "_" gene;
          print $1, col2, "chr", "pos", "variant_id", "pval_nominal"
        }' intdata/102-GTExv10_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt > temp_file.txt && mv temp_file.txt intdata/102-GTExv10_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions_200kb_"${namestring}".txt

done < "intdata/92-regions_of_interest_for_plink_selection.txt"
```

And the eQTL datasets are now assembled using the LocusFocus publicly available python script, into the HTML files that allow for upload to the web tool.

```{bash}
for file in $(ls intdata/102-GTExv10_eQTL_female/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions*.txt); do
    pre=$(basename "$file")
    temp=${pre#*_*_}
    reg=${temp%.*}
    
    out=intdata/94-female_${reg}_snps_eqtlv10.html
    
    temp="${reg/_/:}"
  region="${temp/_/-}"
    
    python3 code/merge_and_convert_to_html.py $file $region $out
done

for file in $(ls intdata/102-GTExv10_eQTL_male/Fixed/subset/Fixed2/split_by_gene/subset2/descriptions*.txt); do
    pre=$(basename "$file")
    temp=${pre#*_*_}
    reg=${temp%.*}
    
    out=intdata/94-male_${reg}_snps_eqtlv10.html
    
    temp="${reg/_/:}"
  region="${temp/_/-}"
    
    python3 code/merge_and_convert_to_html.py $file $region $out
done
```

## Upload to Locus Zoom

For each region (surrounding the 9 peaks highlighted in the manuscript), the sex-specific colocalization test can be done using the Locus Focus web tool, by uploading the trio of files:

* Sex-specific eQTL statistics: intdata/94-${sex}_${reg}_snps_eqtlv10.html
* Sex-specific LD matrix: intdata/94-${sex}_${reg}_snps_ldcorrected.ld
* Sex-specific GWAS statistics: intdata/94-${sex}_${reg}_res.txt