# README

A README for tracking manual analyses using non-outlying individuals determined from PCA. It contains information run on subset of individuals from the original PstI cut site dataset consisting for 24F and 24M wild delta smelt individuals.

## General Setup
```bash
cd /group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-fst
mkdir -p  subset_analyses/maleRef subset_analyses/femaleRef
cp Fst_Results_delta_smelt_both_reference_genomes.R subset_analyses/femaleRef/Fst_results_femaleRef.R
cp Fst_Results_delta_smelt_both_reference_genomes.R subset_analyses/maleRef/Fst_results_maleRef.R
cd subset_analyses
for i in $(echo snpDen_females.txt snpDen_males.txt)
do
    prefix=$(echo ${i} | cut -d. -f1)
    newname=$(echo ${prefix}.subset.txt)
    cp ../../step2-snpdensity/individual_SNP_density/femaleRef/${i} ./femaleRef/${newname}
    cp ../../step2-snpdensity/individual_SNP_density/maleRef/${i} ./maleRef/${newname}
done
```

## Female Reference Genome 

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX47_1_F_G08_ATAGCGAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r

```bash
cd femaleRef/
sed -i '/RX41/d' snpDen_males.subset.txt
sed -i '/RX47/d' snpDen_females.subset.txt
sed -i '/RX59/d' snpDen_males.subset.txt
sed -i '/RX68/d' snpDen_males.subset.txt

# take indiviudals out of vcf file
screen
srun -t 72:00:00 --mem=4Gb -p high --pty /bin/bash
source $HOME/.srun_profile
conda activate bcftools
old_vcf="femaleRef_RAD-Pst1_24F_24M.filtered_PASS.biallelic.vcf.gz"
new_vcf="femaleRef_RAD-Pst1_24F_24M.filtered_PASS.biallelic.subset.vcf.gz"
ln -s ../../../step2-snpdensity/individual_SNP_density/femaleRef/${old_vcf} ${old_vcf}
ln -s ../../../step2-snpdensity/individual_SNP_density/femaleRef/${old_vcf}.tbi ${old_vcf}.tbi
cat snpDen*males.subset.txt >> femaleRef_all_individuals.txt
bcftools view --samples-file femaleRef_all_individuals.txt ${old_vcf} > ${new_vcf}
out_prefix=$(echo ${new_vcf} | cut -d. -f1-4)
vcftools --gzvcf ${new_vcf} --weir-fst-pop snpDen_males.subset.txt --weir-fst-pop snpDen_females.subset.txt --out ${out_prefix}

conda activate difcover # has R w packages loaded
R 
```

See `subset_analyses/maleRef/Fst_results_maleRef.R` for R work...





## Male Reference Genome

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r

```bash


```

