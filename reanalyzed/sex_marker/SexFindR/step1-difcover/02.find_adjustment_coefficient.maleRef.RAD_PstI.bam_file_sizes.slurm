#!/bin/bash -l
#SBATCH -J ACbamM
#SBATCH -e slurm/j%j.find_adjustment_coefficient.maleRef.RAD_PstI.bam_file_sizes.err
#SBATCH -o slurm/j%j.find_adjustment_coefficient.maleRef.RAD_PstI.bam_file_sizes.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G

# Script to make file of BAM file sizes to calculate AC's

###############
###  SETUP  ###
###############


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
tech="RAD"
bam_suffix="sort-n.fixmate-m.sort.markdup-r.bam"
# directories
bam_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-mapping/00-alignments_${sex}"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover/${SEX}Ref/${tech}_stats"
# files
outfile="${SEX}Ref_aligned_bam_file_sizes.list"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    Sex : ${sex}"
echo "    RAD bam Directory : ${bam_dir}" 
echo "    Output Directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# get file sizes
cd ${bam_dir}
ls -l *.${bam_suffix} | awk '{print $9, $5}' > ${out_dir}/${outfile}
# find best of each sex
cd ${out_dir}
female_one=$(sort -n -k 2 -r ${SEX}Ref_aligned_bam_file_sizes.list | grep "_F_" | head -n 1)
female_two=$(sort -n -k 2 -r ${SEX}Ref_aligned_bam_file_sizes.list | grep "_F_" | head -n 2 | tail -n 1)
male_one=$(sort -n -k 2 -r ${SEX}Ref_aligned_bam_file_sizes.list | grep "_M_" | head -n 1)
male_two=$(sort -n -k 2 -r ${SEX}Ref_aligned_bam_file_sizes.list | grep "_M_" | head -n 2 | tail -n 1)
# print to file
echo -e "Top females: \n${female_one}\n${female_two}\n\nTop males: \n${male_one}\n${male_two}" > ${SEX}Ref_aligned_bam_top2.txt
echo -e "Top females: \n${female_one}\n${female_two}\n\nTop males: \n${male_one}\n${male_two}"
echo ""

# get AC's of tops
## make ratios 
ratio_11=$(echo $(echo "scale=4; $(echo ${male_one} | awk '{print $2}')/$(echo ${female_one} | awk '{print $2}')" | bc -l))
ratio_12=$(echo $(echo "scale=4; $(echo ${male_two} | awk '{print $2}')/$(echo ${female_one} | awk '{print $2}')" | bc -l))
ratio_21=$(echo $(echo "scale=4; $(echo ${male_one} | awk '{print $2}')/$(echo ${female_two} | awk '{print $2}')" | bc -l))
ratio_22=$(echo $(echo "scale=4; $(echo ${male_two} | awk '{print $2}')/$(echo ${female_two} | awk '{print $2}')" | bc -l))
## print ratios
echo "" >> ${SEX}Ref_aligned_bam_top2.txt
echo "AC Values (AC=maleFileSize/femaleFileSize):" >> ${SEX}Ref_aligned_bam_top2.txt
echo ${ratio_11} $(echo ${female_one} | awk '{print $1}') $(echo ${male_one} | awk '{print $1}') >> ${SEX}Ref_aligned_bam_top2.txt 
echo ${ratio_12} $(echo ${female_one} | awk '{print $1}') $(echo ${male_two} | awk '{print $1}') >> ${SEX}Ref_aligned_bam_top2.txt
echo ${ratio_21} $(echo ${female_two} | awk '{print $1}') $(echo ${male_one} | awk '{print $1}') >> ${SEX}Ref_aligned_bam_top2.txt
echo ${ratio_22} $(echo ${female_two} | awk '{print $1}') $(echo ${male_two} | awk '{print $1}') >> ${SEX}Ref_aligned_bam_top2.txt
## create AC file
cd ${out_dir}
cat << ACFILE >> ${SEX}Ref_aligned_bam_top2.AC
${ratio_11} $(echo ${female_one} | awk '{print $1}') $(echo ${male_one} | awk '{print $1}')
${ratio_12} $(echo ${female_one} | awk '{print $1}') $(echo ${male_two} | awk '{print $1}')
${ratio_21} $(echo ${female_two} | awk '{print $1}') $(echo ${male_one} | awk '{print $1}')
${ratio_22} $(echo ${female_two} | awk '{print $1}') $(echo ${male_two} | awk '{print $1}')
ACFILE


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
