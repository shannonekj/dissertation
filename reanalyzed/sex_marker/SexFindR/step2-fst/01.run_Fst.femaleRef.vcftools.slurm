#!/bin/bash -l
#SBATCH -J femFst
#SBATCH -e slurm/j%j.run_Fst.femaleRef.vcftools.err
#SBATCH -o slurm/j%j.run_Fst.femaleRef.vcftools.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=4G

# Script to  using female reference genome
# FILTERS:
##    -PASS vcf filter
##    -biallelic sites
# SOFTWARE:
##    -vcftools=0.1.14

###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate bcftools

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"

### MAKE CHANGES ###
# variables
sex="F"
SEX="female"
# directories
snp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-snpdensity/individual_SNP_density/${SEX}Ref"
vcf_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-variantcalling/${SEX}Ref"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-fst"
# files
vcf_file="${SEX}Ref_RAD-Pst1_24F_24M.filtered_PASS.biallelic.vcf.gz"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    SNP Density directory"
echo "    VCF directory : ${vcf_dir}" 
echo "    VCF file : ${vcf_file}"
echo "    output directory : ${out_dir}"
echo ""

# you know the drill
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}
echo $(date +%D' '%T) "Finding Fst between females and males using the  ${SEX} reference genome"
vcftools --gzvcf ${vcf_dir}/${vcf_file} --weir-fst-pop ${snp_dir}/snpDen_males.txt --weir-fst-pop ${snp_dir}/snpDen_females.txt --out biallelic_fst_${SEX}Ref
echo "    done."


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
