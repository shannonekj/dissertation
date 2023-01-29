#!/bin/bash -l
#SBATCH -J mGWAS
#SBATCH -e slurm/j%j.gwas.maleRef.err
#SBATCH -o slurm/j%j.gwas.maleRef.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G

#Script to run gwas using:
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
# directories
vcf_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-variantcalling/${SEX}Ref"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas"
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
echo "    Reference Genome : ${reference}"
echo "    VCF File : $(ls ${vcf_dir}/${vcf_file})"
echo "    GWAS Directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# make metadata file
echo $(date +%D' '%T) "Making metadata file"
[ -f sex_samples_metadata.txt ] && rm sex_samples_metadata.txt
bcftools query -l ${vcf_dir}/${vcf_file} > tmpfile
for line in $(cat tmpfile)
do
    pheno=$(echo ${line} | cut -d_ -f2)
    echo "    ${pheno} ${line}"
    echo -e "${line}\t${line}\t${pheno}" >> sex_samples_metadata.txt
done
rm tmpfile
echo "    done."

# filter VCF & convert to plink fmt
echo $(date +%D' '%T) "Formatting VCF for plink & filtering to remove sites where missingData>50%, maf<5% or maf>95%"
cd ${out_dir}
vcftools --gzvcf ${vcf_dir}/${vcf_file} --plink --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --out ${SEX}Ref_gwas
echo "    done."
echo ""

# run plink 
echo $(date +%D' '%T) "Running plink to generate plink binary files"
plink --file ${SEX}Ref_gwas --pheno sex_samples_metadata.txt --make-bed --out ${SEX}Ref_gwas_plink --noweb --allow-no-sex
echo "    done."
echo ""

# run GEMMA
echo $(date +%D' '%T) "Running GEMMA using a likelihood ratio test"
gemma -bfile ${SEX}Ref_gwas_plink -lm 2 -o ${SEX}Ref_gwas_gemma_out
echo "    done."
echo ""


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
