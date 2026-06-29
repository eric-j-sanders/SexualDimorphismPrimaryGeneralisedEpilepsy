# Purpose

When the codes in this folder are executed, the summary statistics from the sex-stratified GWAS of JME are compared between females and males, to test for shared genetic architecture.

# Input files

* The listed outputs of 05-SexStratifiedGWAS
* Scripts for the GPS-GEV test, copied from https://github.com/twillis209/gps_cpp/tree/master

# Primary outputs

* A list of the SNPs featured in both the male and female GWAS, with Cochrans Q statistics for each: intdata/100-gencor/05-matched_pruned_pvalues.txt
* The printed output of the GPS-GEV test: intdata/100-gencor/07-gps_pvalue.txt
* Manhattan plot of IVW meta-P values for each SNP, colored by whether signal was strongest when sex was considered separately or together: figures/paper/meta_plot.png

# Code

## Collecting P values of SNPs in both GWAS

Now we need to prepare a file that satisfies:

* Only has SNPs present in both studies
* Is LD-pruned to SNPs that aren’t highly correlated
* For each row has p_males p_females

First we will get the list of SNPs present in both study outputs:

In R:

```{r}
library(data.table)

male=fread('intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)
female=fread('intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)

snp_both = intersect(male$ID,female$ID)

fwrite(data.frame(id=snp_both),'intdata/100-gencor/snps_bothgwas.txt',col.names=FALSE)
```

And now we subset the original GWAS data down to only the SNPs that made it through both GWAS, then do LD pruning:

```{bash}
echo "." > intdata/100-gencor/exclude_dot.txt # Specifically excluded any SNPs with ID "."

plink2 \
  --bfile intdata/55-merge_males_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
  --extract intdata/100-gencor/snps_bothgwas.txt \
  --exclude intdata/100-gencor/exclude_dot.txt \
  --max-alleles 2 \
  --make-bed \
  --out intdata/100-gencor/males_gwassnps

plink2 \
  --bfile intdata/55-merge_females_imputed_eur_unrel_matched_maf01both_hwe04_hardcalls \
  --extract intdata/100-gencor/snps_bothgwas.txt \
  --exclude intdata/100-gencor/exclude_dot.txt \
  --max-alleles 2 \
  --make-bed \
  --out intdata/100-gencor/females_gwassnps

echo "intdata/100-gencor/females_gwassnps" > intdata/100-gencor/merge_list.txt

plink \
  --bfile intdata/100-gencor/males_gwassnps \
  --bmerge intdata/100-gencor/females_gwassnps \
  --make-bed \
  --out intdata/100-gencor/merged_males_females_gwassnps

plink2 \
  --bfile intdata/100-gencor/merged_males_females_gwassnps \
  --indep-pairwise 50 5 0.8 \
  --out intdata/100-gencor/04-pruned_snps_for_GPS_GEV
```

In R, now assemble the final data of LD-pruned SNPs for the test, with the male-specific and female-specific p-values:

```{r}
library(data.table)

keep_snps = fread("intdata/100-gencor/04-pruned_snps_for_GPS_GEV.prune.in", header = FALSE)$V1

male=fread('intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)
female=fread('intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt',data.table=FALSE)

male_matched = male[male$ID %in% keep_snps,]
male_matched = male_matched[order(male_matched$ID),]
female_matched = female[female$ID %in% keep_snps,]
female_matched = female_matched[order(female_matched$ID),]

identical(male_matched$ID,female_matched$ID)

matched_all = data.frame(rsid = male_matched$ID,
                  p_male = male_matched$P,
                  or_male = male_matched$OR,
                  se_male = abs(log(male_matched$OR)) / qnorm(male_matched$P/2, lower.tail = FALSE),
                  p_female = female_matched$P,
                  or_female = female_matched$OR,
                  se_female = abs(log(female_matched$OR)) / qnorm(female_matched$P/2, lower.tail = FALSE))

matched_all = matched_all[!is.na(matched_all$se_male) & !is.na(matched_all$se_female), ]

matched_all = matched_all[
  !is.na(matched_all$se_male) & 
  !is.na(matched_all$se_female) &
  matched_all$se_male > 0 &
  matched_all$se_female > 0, ]
```

## Calculate Q Statistics and Heterogeneity Lambda GC

In R, continue from last code block, preserving its defined object environment:

```{r}
w_male = 1 / matched_all$se_male^2
w_female = 1 / matched_all$se_female^2

pooled = (w_male * log(matched_all$or_male) + w_female * log(matched_all$or_female)) / 
          (w_male + w_female)

matched_all$Q = w_male   * (log(matched_all$or_male)   - pooled)^2 +
                 w_female * (log(matched_all$or_female) - pooled)^2

matched_all$Q_pval = pchisq(matched_all$Q, df = 1, lower.tail = FALSE)

summary(matched_all$Q)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.00000  0.09924  0.44463  0.97128  1.29086 24.91625
summary(matched_all$Q_pval)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000006 0.2558894 0.5048976 0.5038589 0.7527432 1.0000000 

median(qchisq(1-matched_all$Q_pval,2))/qchisq(0.5,2) # lambda gc = 0.9859372

fwrite(matched_all,
       "intdata/100-gencor/05-matched_pruned_pvalues.txt",sep='\t')
```

## GPS-GEV Test

In Bash:

```{bash}
~/Programs/gps_cpp/build/apps/computeGpsCLI -i intdata/100-gencor/05-matched_pruned_pvalues.txt -o intdata/100-gencor/06-Gps_output.txt -a p_male -b p_female -c Male -d Female

sbatch --mem=120G --cpus-per-task=4 --time=20:00:00 --wrap="~/Programs/gps_cpp/build/apps/permuteTraitsCLI -i intdata/100-gencor/05-matched_pruned_pvalues.txt -o intdata/100-gencor/06-Gps_nullsim_output.txt -a p_male -b p_female -d 3000"

# Wait for null permuted test statistics to be collected, then proceed:

Rscript ~/Programs/gps_cpp/R/fit_gevd_and_compute_pvalue.R \
  -g 20.2837 \
  -p intdata/100-gencor/06-Gps_nullsim_output.txt \
  -o intdata/100-gencor/07-gps_pvalue.txt
```

## Inverse-Variance Weighted Meta-Analysis

In R:

```{r}
library(dplyr)
library(data.table)
library(ggplot2)

male = fread('intdata/91-males_forLZ_maf_hwe_sfscf_pipeline.txt',select  = c("#CHROM","POS","ID","REF","ALT", "A1", "OR", "P"),data.table=FALSE)
female = fread('intdata/91-females_forLZ_maf_hwe_sfscf_pipeline.txt',select  = c("#CHROM","POS","ID","REF","ALT", "A1", "OR", "P"),data.table=FALSE)

male$A2 = ifelse(male$A1 == male$ALT, male$REF, male$ALT)
female$A2 = ifelse(female$A1 == female$ALT, female$REF, female$ALT)

male$REF = NULL
male$ALT = NULL
female$REF = NULL
female$ALT = NULL

both = inner_join(male,female,by=c('#CHROM','POS','ID','A1','A2'),suffix=c('.male','.female'))

both$beta.male = log(both$OR.male)
both$beta.female = log(both$OR.female)

both$se.male = abs(both$beta.male/qnorm(both$P.male/2))
both$se.female = abs(both$beta.female/qnorm(both$P.female/2))

both$beta.meta = (both$beta.male/both$se.male^2 + both$beta.female/both$se.female^2)/(1/both$se.male^2 + 1/both$se.female^2)

both$se.meta = sqrt(1/(1/both$se.male^2 + 1/both$se.female^2))

both$p.meta = 2*pnorm(-abs(both$beta.meta/both$se.meta))

both$lowest = ifelse(both$P.male < both$P.female & both$P.male < both$p.meta, "male",
               ifelse(both$P.female < both$p.meta, "female", "meta"))

both$minp = pmin(both$P.male,both$P.female,both$p.meta)

both01 = both[which(both$p.meta<0.1),]

both01$`#CHROM` = factor(both01$`#CHROM`,levels=c(1:22,"X"))

chr_offsets = both01 %>%
  group_by(`#CHROM`) %>%
  summarise(chr_len = max(POS)) %>%
  mutate(offset = cumsum(as.numeric(lag(chr_len, default = 0))))

both01 = both01 %>%
  left_join(chr_offsets, by = "#CHROM") %>%
  mutate(pos_cum = as.numeric(POS) + offset)

axis_df = both01 %>%
  group_by(`#CHROM`) %>%
  summarise(center = (max(pos_cum) + min(pos_cum)) / 2)

chr_boundaries = both01 %>%
  group_by(`#CHROM`) %>%
  summarise(boundary = min(pos_cum)) %>%
  filter(as.integer(`#CHROM`) > 1)

both01 = both01[sample(nrow(both01)),]

snps_to_label = c('rs17318744','rs191442248','rs702041','rs806795','rs4938014','rs8058223','rs4889844','rs2525035','rs197013')

label_df = both01[both01$ID %in% snps_to_label, ]



mygg = ggplot(both01, aes(x = pos_cum, y = -log10(p.meta), color = lowest))  +
  geom_point(size = 1, alpha = 0.6)+
  geom_vline(data = chr_boundaries, aes(xintercept = boundary), linewidth = 0.3) +
  scale_color_manual(
    values = c("male" = "#2166ac", "female" = "#d6604d", "meta" = "#4dac26"),
    name = "Lowest P"
  ) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$`#CHROM`) +
  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "black") +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_bw(base_size=16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

mygg = mygg +
  geom_point(data = label_df, aes(x = pos_cum, y = -log10(p.meta)),color='black', size = 4)+
  geom_point(data = label_df, aes(x = pos_cum, y = -log10(p.meta),color=lowest), size = 3) + 
  geom_label(data = label_df, aes(x = pos_cum, y = -log10(p.meta), label = ID),
                   color = "black", size = 6,nudge_x=1.8e8)

ggsave('figures/paper/meta_plot.png',mygg,width=18,height=9,units='in',dpi=300)
```