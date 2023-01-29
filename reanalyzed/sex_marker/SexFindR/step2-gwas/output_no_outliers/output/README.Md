# README for manual processing

# Step 1: Pull data into R
bash
```
conda activate gwas
cd /group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas/output_no_outliers/output
R
```

R
```
setwd("/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas/output_no_outliers/output")
library(tidyverse)

## female
fdat <- read_tsv("femaleRef_no_outliers_gemma_out.assoc.txt", col_names=T)
fdat_sub_splt <- fdat %>% filter(p_lrt <= 0.005) %>% separate(rs, c("chromosome", "position"), sep=":", remove=FALSE)
write_tsv(fdat_sub_splt, "femaleRef_no_outliers_gemma_out.pval5e-3.assoc.txt")

## male
mdat <- read_tsv("maleRef_no_outliers_gemma_out.assoc.txt", col_names=T)
mdat_sub_splt <- mdat %>% filter(p_lrt <= 0.005) %>% separate(rs, c("chromosome", "position"), sep=":", remove=FALSE)
write_tsv(mdat_sub_splt, "maleRef_no_outliers_gemma_out.pval5e-3.assoc.txt")
```

bash
```
tail -n +2 femaleRef_no_outliers_gemma_out.pval5e-3.assoc.txt | awk '{print $2}' >> femaleRef_no_outliers_gemma_out.pval5e-3.loci
tail -n +2 maleRef_no_outliers_gemma_out.pval5e-3.assoc.txt | awk '{print $2}' >> maleRef_no_outliers_gemma_out.pval5e-3.loci
conda activate gwas
plink --bfile ../femaleRef_no_outliers_gemma_in --extract femaleRef_no_outliers_gemma_out.pval5e-3.loci --recode compound-genotypes --allow-extra-chr --out femaleRef_no_outliers_gemma_out.genotypes.pval5e-3
plink --bfile ../maleRef_no_outliers_gemma_in --extract maleRef_no_outliers_gemma_out.pval5e-3.loci --recode compound-genotypes --allow-extra-chr --out maleRef_no_outliers_gemma_out.genotypes.pval5e-3
```
