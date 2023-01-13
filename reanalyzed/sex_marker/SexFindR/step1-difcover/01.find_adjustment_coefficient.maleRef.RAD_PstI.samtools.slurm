#!/bin/bash -l
#SBATCH -J ACbamM
#SBATCH -e slurm/j%j.find_adjustment_coefficient.maleRef.RAD_PstI.samtools.err
#SBATCH -o slurm/j%j.find_adjustment_coefficient.maleRef.RAD_PstI.samtools.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=16G

# Script to find adjustment coefficient between two samples' modal coverage 
# SOFTWARE:
##    -samtools=1.10

###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate samtools

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SOFTWARE VERSIONS:"
samtools --version

### MAKE CHANGES ###
sex="M"
SEX="male"
tech="RAD"
bam_suffix="sort-n.fixmate-m.sort.markdup-r.bam"
# directories
bam_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-mapping/00-alignments_${sex}"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover/${SEX}Ref/${tech}_stats"
# files
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
echo "    RAD bam Directory : ${bam_dir}" 
echo "    Output Directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# loop through to get per sample coverage distribution
for i in ${bam_dir}/*.${bam_suffix}
do
echo "$(date +%D' '%T) Processing ${i}"
filename=$(echo ${i} | rev | cut -d/ -f1 | rev)
echo "    getting stats... $(date +%T)"
samtools stats ${i} > ${filename}.stats
echo "    adding stat to master list"
number=$(grep ^COV ${filename}.stats | cut -f 2- | awk -v max=0 '{if($3>max){want=$2; max=$3}}END{print want} ' -)
echo $filename $number >> ALL_modal_depth.txt
echo "    done."
done

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
