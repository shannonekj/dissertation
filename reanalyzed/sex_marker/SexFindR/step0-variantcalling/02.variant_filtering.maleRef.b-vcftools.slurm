#!/bin/bash -l
#SBATCH -J fltMvcf
#SBATCH -e slurm/j%j.variant_filtering.maleRef.platypus.err
#SBATCH -o slurm/j%j.variant_filtering.maleRef.platypus.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G

# Script to filter variants called from RAD reads aligned w/samtools using male reference genome
# FILTERS:
##    -PASS vcf filter
##    -biallelic sites
# SOFTWARE:
##    -bcftools=1.9
##    -vcftools=0.1.14

###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate platypus # this has [b/v]cftools installed 

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SOFTWARE VERSIONS:"
bcftools --version
vcftools --version

### MAKE CHANGES ###
sex="M"
SEX="male"
# directories
bam_dir="/home/sejoslin/git_repos/SexFindR/analysis_2022/step0-mapping/00-alignments_${sex}"
out_dir="/home/sejoslin/git_repos/SexFindR/analysis_2022/step0-variantcalling/${SEX}Ref"
# files
vcf_prefix="${SEX}Ref_RAD-Pst1_24F_24M"
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
echo "    Reference Genome : ${reference}"
echo "    vcf file : $(ls ${out_dir}/${vcf_prefix}.vcf)"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# filter for PASS
echo $(date +%D' '%T) "filter for variants that PASS"
bcftools view -f PASS ${vcf_prefix}.vcf > ${out_dir}/${vcf_prefix}.filtered_PASS.vcf
echo "    done."
echo ""

# filter for biallelic sites
echo $(date +%D' '%T) "filter for biallelic sites"
cd ${out_dir}
vcftools --vcf ${vcf_prefix}.filtered_PASS.vcf --max-alleles 2 --stdout --recode --recode-INFO-all | gzip -c > ${vcf_prefix}.filtered_PASS.biallelic.vcf.gz
echo "    done."
echo ""

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
