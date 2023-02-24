# README

A README for tracking manual analyses using non-outlying individuals determined from PCA. It contains information run on subset of individuals from the original PstI cut site dataset consisting for 24F and 24M wild delta smelt individuals.

Filters:
- rm outliers
- rm pruned sites from GWAS
- rm loci missing>50%, 5%<maf<95%, indels

## General Setup
```bash
cd /group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-fst
mkdir -p  subset_analyses/maleRef subset_analyses/femaleRef
cp Fst_Results_delta_smelt_both_reference_genomes.R subset_analyses/Fst_results.subset_analysis.R
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

# make vcf file from GWAS plink files
screen
srun -t 72:00:00 --mem=4Gb -p high --pty /bin/bash
source $HOME/.srun_profile
conda activate gwas
sex="female"
plink_dir="../../../step2-gwas/output_no_outliers"
plink_prefix="${sex}Ref_no_outliers_gemma_in"
for i in $(echo fam bim bed)
do
    ln -s ${plink_dir}/${plink_prefix}.${i} ${plink_prefix}.${i}
done
plink --bfile ${plink_prefix} --double-id --allow-no-sex --allow-extra-chr --recode vcf-iid --out ${sex}Ref_no_outliers_pruned

# index vcf 
bgzip -c ${sex}Ref_no_outliers_pruned.vcf > ${sex}Ref_no_outliers_pruned.vcf.gz
tabix -p vcf ${sex}Ref_no_outliers_pruned.vcf.gz

# run fst
vcftools --gzvcf ${sex}Ref_no_outliers_pruned.vcf.gz --weir-fst-pop snpDen_males.subset.txt --weir-fst-pop snpDen_females.subset.txt --out ${sex}Ref_no_outliers_pruned


## now run maleRef analysis

##OLD
#conda activate difcover # has R w packages loaded
#R 
```

#See `subset_analyses/femaleRef/Fst_results_femaleRef.R` for R work...





## Male Reference Genome

Excluded individuals:
* RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
* RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
* RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r

```bash
cd /group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-fst/subset_analyses/maleRef
sed -i '/RX41/d' snpDen_males.subset.txt
sed -i '/RX59/d' snpDen_males.subset.txt
sed -i '/RX68/d' snpDen_males.subset.txt

sex="male"
plink_dir="../../../step2-gwas/output_no_outliers"
plink_prefix="${sex}Ref_no_outliers_gemma_in"
for i in $(echo fam bim bed)
do
    ln -s ${plink_dir}/${plink_prefix}.${i} ${plink_prefix}.${i}
done
plink --bfile ${plink_prefix} --double-id --allow-no-sex --allow-extra-chr --recode vcf-iid --out ${sex}Ref_no_outliers_pruned

# index vcf
bgzip -c ${sex}Ref_no_outliers_pruned.vcf > ${sex}Ref_no_outliers_pruned.vcf.gz
tabix -p vcf ${sex}Ref_no_outliers_pruned.vcf.gz

# run fst
vcftools --gzvcf ${sex}Ref_no_outliers_pruned.vcf.gz --weir-fst-pop snpDen_males.subset.txt --weir-fst-pop snpDen_females.subset.txt --out ${sex}Ref_no_outliers_pruned

cd ../
conda activate difcover
R
```

See `Fst_results.subset_analyses.R` for R work




