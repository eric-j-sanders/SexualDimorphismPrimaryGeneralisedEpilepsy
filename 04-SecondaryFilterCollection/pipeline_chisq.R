library(data.table)
library(dplyr)
library(qqman)
library(fuzzyjoin)

for_removal=data.frame()





# CASES

in.batch = 'intdata/86-freq_batchimpute_cases.afreq'
in.joint = 'intdata/86-freq_jointimpute_cases.afreq'
out.diff = 'intdata/86-result_testing_af_by_pipeline_CASES.txt'

freq.b = fread(in.batch,data.table=FALSE)
freq.j = fread(in.joint,data.table=FALSE)

freq.j$`#CHROM` = as.character(freq.j$`#CHROM`)

# Merge and only keep SNPs in both imputation files
res = full_join(freq.b[,c(1:5,7,8)], freq.j[,c(1:5,7,8)], by = c('#CHROM','POS','REF','ALT'), suffix=c('_batch','_joint'))

res.both = res[complete.cases(res),]

# Start to assemble as allele counts for testing
res.both$ALT_CT_batch = round(res.both$OBS_CT_batch*res.both$ALT_FREQS_batch)
res.both$ALT_CT_joint = round(res.both$OBS_CT_joint*res.both$ALT_FREQS_joint)


compute_chisq_p = function(alt1,total1,alt2,total2){
	ref1 = total1-alt1
	ref2 = total2-alt2
	
	tab = matrix(c(alt1, ref1, alt2, ref2), nrow = 2, byrow = TRUE)
	
	return(chisq.test(tab, correct = FALSE)$p.value)
}

res.both$chisq_p = mapply(compute_chisq_p, 
													res.both$ALT_CT_batch, 
													res.both$OBS_CT_batch, 
													res.both$ALT_CT_joint, 
													res.both$OBS_CT_joint)

# Save the case-specific SNPs to remove
fwrite(res.both,file=out.diff,sep='\t')

# Keep a data frame of these, to add to later when case and control lists are merged:
for_removal = rbind(for_removal,res.both)
for_removal$group = "CASES"

# Tidy up environment before moving on to controls
rm(list = setdiff(ls(), "for_removal"))
gc()

# CONTROLS
# (Following code basically copies all the previous code)

in.batch = 'intdata/86-freq_batchimpute_controls.afreq'
in.joint = 'intdata/86-freq_jointimpute_controls.afreq'
out.diff = 'intdata/86-result_testing_af_by_pipeline_CONTROLS.txt'

freq.b = fread(in.batch,data.table=FALSE)
freq.j = fread(in.joint,data.table=FALSE)

freq.j$`#CHROM` = as.character(freq.j$`#CHROM`)

res = full_join(freq.b[,c(1:5,7,8)], freq.j[,c(1:5,7,8)], by = c('#CHROM','POS','REF','ALT'), suffix=c('_batch','_joint'))

res.both = res[complete.cases(res),]

res.both$ALT_CT_batch = round(res.both$OBS_CT_batch*res.both$ALT_FREQS_batch)
res.both$ALT_CT_joint = round(res.both$OBS_CT_joint*res.both$ALT_FREQS_joint)

compute_chisq_p = function(alt1,total1,alt2,total2){
	ref1 = total1-alt1
	ref2 = total2-alt2
	
	tab = matrix(c(alt1, ref1, alt2, ref2), nrow = 2, byrow = TRUE)
	
	return(chisq.test(tab, correct = FALSE)$p.value)
}
res.both$chisq_p = mapply(compute_chisq_p, 
													res.both$ALT_CT_batch, 
													res.both$OBS_CT_batch, 
													res.both$ALT_CT_joint, 
													res.both$OBS_CT_joint)

fwrite(res.both,file=out.diff,sep='\t')

res.both$group = "CONTROLS"

for_removal = rbind(for_removal,res.both)

# Save the data set that has SNPs differing in EITHER cases or controls:
fwrite(for_removal,file='intdata/86-variants_with_differing_af_by_pipeline_CASES_OR_CONTROLS.txt')