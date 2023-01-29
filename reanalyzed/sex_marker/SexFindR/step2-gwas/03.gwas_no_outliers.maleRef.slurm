#!/bin/bash -l
#SBATCH -J mGWAS
#SBATCH -e slurm/j%j.gwas_no_outliers.maleRef.err
#SBATCH -o slurm/j%j.gwas_no_outliers.maleRef.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G


# Script to run gwas on individuals passing filtration:
# REMOVED INDIVIDUALS:
##    -RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
##    -RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
##    -RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r


# SOFTWARE:
##    -bcftools=1.9
##    -htslib=1.9
##    -plink=v1.90b6.21
##    -vcftools=0.1.14
##    -GEMMA=0.98.3

###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh
conda activate gwas

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"

### MAKE CHANGES ###
sex="M"
SEX="male"
prefix="${SEX}Ref_no_outliers"
# directories
vcf_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-variantcalling/${SEX}Ref"
gwa_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas"
out_dir="${gwa_dir}/output_no_outliers"
# files
vcf_file="${SEX}Ref_RAD-Pst1_24F_24M.filtered_PASS.biallelic.vcf.gz"
if [ $sex == "M" ]
then
    reference="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/04-linkage_map/out_chromonomer/CHRR_integrated.fa"
elif [ $sex == "F" ]
then
    reference="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/03-assemblies/sandbox_ipa_F/Hyp_tra_F_combined_chromonomer/output/CHRR_integrated.fa"
else
    echo "sex input is not female or male. exiting process."
    exit
fi
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    Sex : ${sex}"
echo "    Prefix : ${prefix}"
echo "    Reference Genome : ${reference}"
echo "    VCF File : $(ls ${vcf_dir}/${vcf_file})"
echo "    GWAS Directory : ${gwa_dir}"
echo "    Output Directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# make outlier list
[ -f ${SEX}Ref_outliers.list ] && rm ${SEX}Ref_outliers.list
cat << EOT >> ${SEX}Ref_outliers.list
RX41_2_M_F11_AAGGACAC.sort-n.fixmate-m.sort.markdup-r
RX59_2_M_H11_GATGAATC.sort-n.fixmate-m.sort.markdup-r
RX68_2_M_G12_CACCTTAC.sort-n.fixmate-m.sort.markdup-r
EOT
## START HERE

# plink things
VCF="${vcf_dir}/${vcf_file}"
## filter VCF
echo $(date +%D' '%T) "Formatting VCF for plink & filtering to remove sites where missingData>50%, maf<5% or maf>95%"
cd ${out_dir}
call="vcftools --gzvcf ${VCF} --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --recode --recode-INFO-all --stdout | gzip -c > ${prefix}.vcf.gz"
echo ${call}
eval ${call}
echo "    done."
echo ""

## linkage prune
echo $(date +%D' '%T) "Calculate R2 & identify linked SNPs (R2>0.1)"
call="plink --vcf ${prefix}.vcf.gz --double-id --remove-fam ${SEX}Ref_outliers.list --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ${prefix}"
echo ${call}
eval ${call}
echo "    done."
echo ""

## run pca just to be sure
### this also makes BED
echo $(date +%D' '%T) "Perform PCA & make bed"
call="plink --vcf ${prefix}.vcf.gz --double-id --remove-fam ${SEX}Ref_outliers.list --allow-extra-chr --set-missing-var-ids @:# --extract ${prefix}.prune.in --allow-no-sex --make-bed --pca --out ${prefix}"
echo ${call}
eval ${call}
echo "    done."
echo ""



# GEMMA
## make phenotype file
echo $(date +%D' '%T) "Creating phenotype file"
[ -f ${SEX}Ref_phenotypes.list ] && rm ${SEX}Ref_phenotypes.list
bcftools query -l ${VCF} > tmpfile
for line in $(cat tmpfile)
do
    pheno=$(echo ${line} | cut -d_ -f2)
    echo "    ${pheno} ${line}"
    echo -e "${line}\t${line}\t${pheno}" >> ${SEX}Ref_phenotypes.list
done
rm tmpfile
echo "    done."
echo ""

## convert BED to PED + MAP
echo $(date +%D' '%T) "Converting bed file to ped and map"
call="plink --bfile ${prefix} --pheno ${SEX}Ref_phenotypes.list --allow-no-sex --allow-extra-chr --recode tab --out ${prefix}_gemma_in"
echo ${call}
eval ${call}
echo "    done."
echo ""

## make BED, FAM, BIM files with phenotypes
echo $(date +%D' '%T) "Making bed, fam, bim w/ phenotype data"
call="plink --file ${prefix}_gemma_in --pheno ${SEX}Ref_phenotypes.list --remove-fam ${SEX}Ref_outliers.list --extract ${prefix}.prune.in --make-bed --allow-no-sex --allow-extra-chr --out ${prefix}_gemma_in"
echo ${call}
eval ${call}
echo "    done."
echo ""

## likelihood test
echo $(date +%D' '%T) "Running GEMMA using a likelihood ratio test"
call="gemma -bfile ${prefix}_gemma_in -lm 2 -o ${prefix}_gemma_out"
echo ${call}
eval ${call}
echo "    done."
echo ""


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
