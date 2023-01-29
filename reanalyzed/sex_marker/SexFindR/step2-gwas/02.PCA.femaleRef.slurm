#!/bin/bash -l
#SBATCH -J fPCA
#SBATCH -e slurm/j%j.PCA.femaleRef.err
#SBATCH -o slurm/j%j.PCA.femaleRef.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G

#Script to run PCA using:
# SOFTWARE:
##    -plink=v1.90b6.21

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
sex="F"
SEX="female"
prefix="${SEX}Ref_pca"
# directories
vcf_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-variantcalling/${SEX}Ref"
gwa_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas"
out_dir="${gwa_dir}/pca"
# files
vcf_file="${SEX}Ref_RAD-Pst1_24F_24M.filtered_PASS.biallelic.vcf.gz"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    Sex : ${sex}"
echo "    VCF File : $(ls ${vcf_dir}/${vcf_file})"
echo "    GWAS Directory : ${gwa_dir}"
echo "    PCA Directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# plink things
## linkage prune 
VCF="${vcf_dir}/${vcf_file}"
echo $(date +%D' '%T) "Calculate R2 & identify linked SNPs (R2>0.1)"
call="plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ${prefix}"
echo ${call}
eval ${call}
echo "    done."
echo ""

## run PCA
echo $(date +%D' '%T) "Perform PCA"
call="plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${prefix}.prune.in --make-bed --pca --out ${prefix}"
echo ${call}
eval ${call}
echo "    done."
echo ""


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
