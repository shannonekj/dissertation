# README

This directory is an analysis of a subset of individuals from the original PstI cut site dataset consisting for 24F and 24M wild delta smelt individuals.

## General Setup

bash
```
cd /group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-snpdensity/
mkdir -p pruned_dataset/femaleRef/snpden pruned_dataset/maleRef/snpden
cp individual_SNP_density/femaleRef/manual_analysis/SNPdensity_femaleRef_sex.R pruned_dataset/femaleRef/.
cp individual_SNP_density/maleRef/manual_analysis/SNPdensity_maleRef_sex.R pruned_dataset/maleRef/.
```


## Female Reference Genome

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX47_1_F_G08_ATAGCGAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r

bash
```
cd femaleRef/snpden
for i in ../../../individual_SNP_density/femaleRef/*.snpden; do filename=$(echo $i | rev | cut -d/ -f1 | rev); ln -s $i $filename; done
rm Male_RX41_2_M_F11_AAGGACAC.vcf.snpden Female_RX47_1_F_G08_ATAGCGAC.vcf.snpden Male_RX59_2_M_H11_GATGAATC.vcf.snpden Male_RX68_2_M_G12_CACCTTAC.vcf.snpden
cd ../

screen
srun -t 48:00:00 --mem=4Gb -p bigmemh --pty /bin/bash
source $HOME/.srun_profile
conda activate difcover
R
```

For R work see `SNPdensity_femaleRef_sex.R`





## Male Reference Genome

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r
