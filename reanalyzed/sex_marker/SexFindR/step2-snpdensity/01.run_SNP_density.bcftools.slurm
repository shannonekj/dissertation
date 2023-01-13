#!/bin/bash -l
#SBATCH -J cntKmers
#SBATCH -e slurm/j%j.run_SNP_density.bcftools.err
#SBATCH -o slurm/j%j.run_SNP_density.bcftools.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=4G


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
# directories
stp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-snpdensity"
vcf_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-variantcalling"
out_dir="${stp_dir}/individual_SNP_density"
# files
vcf_file="RAD-Pst1_24F_24M.filtered_PASS.biallelic.vcf"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    VCF directory : ${vcf_dir}" 
echo "    VCF file : ${vcf_file}"
echo "    output directory : ${out_dir}"
echo ""

# create run script
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}
[ -f get_SNPdensity.slurm ] && rm get_SNPdensity.slurm
cat << EOT >> get_SNPdensity.slurm
#!/bin/bash
#
#SBATCH -J indVCF
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-02:59:59
#SBATCH --mem=4Gb
#SBATCH -o slurm/j%j.get_SNPdensity.out
#SBATCH -e slurm/j%j.get_SNPdensity.err

# \$1 is name of file containing sample list
# \$2 is output vcf name
# \$3 is input vcf name
# \$4 is sex

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate bcftools

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

echo SCRIPT FOR : \$1

echo "Step 1 : creating individual vcf file"
bcftools view -a -s \$1 -o individual_vcfs/\$4_\$2 \$3

echo "Step 2 : compressing"
bgzip -c individual_vcfs/\$4_\$2 > individual_vcfs/\$4_\${2}.gz

echo "Step 3 : indexing"
bcftools index individual_vcfs/\$4_\${2}.gz

echo "Step 4 : getting SNP density"
vcftools --vcf individual_vcfs/\$4_\$2 --SNPdensity 10000 --out \$4_\$2
EOT
chmod a+x get_SNPdensity.slurm


# make individual vcf files
for i in $(echo -e "femaleRef\nmaleRef")
do
    cd ${out_dir}
    echo $(date +%D' '%T) "Processing ${i}"
    [ -d ${i} ] || mkdir -p ${i}
    cd ${i}
    [ -f get_SNPdensity.slurm ] || cp ../get_SNPdensity.slurm .
    [ -d slurm ] || mkdir slurm
    [ -f ${i}_${vcf_file}.gz ] || cp ${vcf_dir}/${i}/${i}_${vcf_file}.gz ${i}_${vcf_file}.gz
    echo $(date +%D' '%T) "Indexing ${i} vcf file."
    [ -f ${i}_${vcf_file}.gz.tbi ] || gunzip ${i}_${vcf_file}.gz
    [ -f ${i}_${vcf_file}.gz.tbi ] || bgzip -c ${i}_${vcf_file} > ${i}_${vcf_file}.gz
    [ -f ${i}_${vcf_file}.gz.tbi ] || bcftools index -t ${i}_${vcf_file}.gz
    echo "    done."
    mkdir -p individual_vcfs
    for sample in `bcftools query -l ${i}_${vcf_file}.gz`
    do
        echo $(date +%D' '%T) "Making files of individuals for each sex"
        sex=$(echo ${sample} | cut -d_ -f3)
        if [ $sex == "F" ]
        then
            echo ${sample} >> snpDen_females.txt
            SEX="Female"
        elif [ $sex == "M" ]
        then
            echo ${sample} >> snpDen_males.txt
            SEX="Male"
        else
            echo "sex is not male or female! exiting."
            exit
        fi
        sbatch get_SNPdensity.slurm ${sample} ${sample/.sort-n.fixmate-m.sort.markdup-r*/.vcf} ${i}_${vcf_file}.gz ${SEX}
    done
    cd ${vcf_dir}
done


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
