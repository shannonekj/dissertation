#!/bin/bash -l
#
#SBATCH --array=1-4
#SBATCH -J pltFdc
#SBATCH -e slurm/j%j.plot_results.femaleRef.RAD_PstI.R.a%a.err
#SBATCH -o slurm/j%j.plot_results.femaleRef.RAD_PstI.R.a%a.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=03:04:05
#SBATCH --mem=4G

# Script to run DifCover on RAD BAM files

###############
###  SETUP  ###
###############
# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate difcover

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"


### VARIABLES ###

sex="F"
SEX="female"
tech="RAD"
array=${SLURM_ARRAY_TASK_ID}

## difcover 
ac=""
a=10        # minimum coverage for sample1
A=219       # maximum coverage for sample1
b=10        # minimum coverage for sample2
B=240       # maximum coverage for sample2
v=1000      # target number of valid bases in the window
l=500       # minimum size of window to output
p=0.7369656     # enrichment scores threshold (when p=2 will report regions with coverage in sample1 being roughly 4 times larger than coverage in sample2)

# directories
dif_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover"
sex_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover/${SEX}Ref"
sts_dir="${sex_dir}/${tech}_stats"
wrk_dir="${sex_dir}/${tech}_analysis/a${array}"
out_dir="${wrk_dir}/outputs"

[ -d ${sts_dir} ] || echo "WARNING: There is no stats directory...backtrack and get some stats."
cd ${sts_dir}
ac_file="${SEX}Ref_aligned_bam_top2.AC"
ac=$(sed -n ${array}p ${ac_file} | awk '{print $1}')

# files
dnacopyout="${wrk_dir}/sample1_sample2.ratio_per_w_CC0_a${a}_A${A}_b${b}_B${B}_v${v}_l${l}.log2adj_${ac}.DNAcopyout"
## genome
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
### symlink 
ref_short=$(echo $reference | rev | cut -f1 -d/ | rev)
ref_pre=$(echo $ref_short | cut -f1 -d.)_${sex}.fa
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}
[ -h ${ref_pre} ] || ln -s ${reference} ${ref_pre}
scaffold_length_file="${ref_short}_${sex}.scaffold_lengths.txt"
###  END  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    Sex : ${sex}"
echo "    Genome : ${ref_pre}"
echo "    Output Directory : ${out_dir}"
echo "    Stats Directory : ${sts_dir}"
echo ""

# make genome scaffold length file
## FMT : [ scaffold | length ] (no header)
echo "$(date +%D' '%T) Checking for scaffold length file & creating if not present"
[ -f ${scaffold_length_file} ] || cat ${ref_pre} | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > ${scaffold_length_file}
echo "   done."

# run R script to make plots
echo "$(date +%D' '%T) Plotting"
call="Rscript ${dif_dir}/plot_difcover.R ${dnacopyout} ${scaffold_length_file} ${out_dir} ${p} ${array}"
echo "    CMD: " ${call}
eval ${call}
echo "    done."

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
