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
cd femaleRef
ln -s /group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/03-assemblies/sandbox_ipa_F/Hyp_tra_F_combined_chromonomer/output/CHRR_integrated.fa CHRR_integrated_F.fa
cd snpden
for i in ../../../individual_SNP_density/femaleRef/*.snpden; do filename=$(echo $i | rev | cut -d/ -f1 | rev); ln -s $i $filename; done
rm Male_RX41_2_M_F11_AAGGACAC.vcf.snpden Female_RX47_1_F_G08_ATAGCGAC.vcf.snpden Male_RX59_2_M_H11_GATGAATC.vcf.snpden Male_RX68_2_M_G12_CACCTTAC.vcf.snpden
cd ../

screen
srun -t 48:00:00 --mem=4Gb -p bigmemh --pty /bin/bash
source $HOME/.srun_profile
conda activate difcover
samtools faidx CHRR_integrated_F.fa
R
```

For R work see `SNPdensity_femaleRef_sex.R`


## Male Reference Genome

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r

bash
```
cd maleRef
ln -s /group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/04-linkage_map/out_chromonomer/CHRR_integrated.fa CHRR_integrated_M.fa
cd snpden
for i in ../../../individual_SNP_density/maleRef/*.snpden; do filename=$(echo $i | rev | cut -d/ -f1 | rev); ln -s $i $filename; done
rm Male_RX41_2_M_F11_AAGGACAC.vcf.snpden Male_RX59_2_M_H11_GATGAATC.vcf.snpden Male_RX68_2_M_G12_CACCTTAC.vcf.snpden
cd ../

screen
srun -t 48:00:00 --mem=4Gb -p bigmemh --pty /bin/bash
source $HOME/.srun_profile
conda activate difcover
samtools faidx CHRR_integrated_M.fa
R
```

For R work see SNPdensity_maleRef_sex.R
