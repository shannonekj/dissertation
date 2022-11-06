#!/bin/bash
#
#SBATCH -J fstqcSx
#SBATCH -p bigmemh
#SBATCH -t 2-12:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -e slurm/fastqc_RAD.j%j.err
#SBATCH -o slurm/fastqc_RAD.j%j.out

# FASTQC RAD PstI sex marker reads

###############
###  SETUP  ###
###############
# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh
# activate conda env
conda activate fastqc
# fail on weird errors
set -o nounset
set -o errexit

# variables
## independent
THREADS=${SLURM_NTASKS}

# directories
work_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/raw_stats"
data_dir="/group/millermrgrp2/shannon/projects/DS_sex-marker/00-RAD-data"
out_dir="${work_dir}/fastQC"


###################
###  DO THINGS  ###
###################
# fastQC R1 & R2 files
echo $(date +%D' '%T) "Checking to make sure ${out_dir} exists"
[ -d ${out_dir} ] || mkdir -p ${out_dir}
echo $(date +%D' '%T) "    ...it does."
cd ${out_dir}

echo $(date +%D' '%T) "FastQCing thangzZ"
fastqc -t ${THREADS} -f fastq -o ${out_dir} ${data_dir}/*.fastq
echo $(date +%D' '%T) "    Done!"
echo ""

# print job information
env | grep SLURM                      # Print values of current jobs SLURM environment variables
scontrol show job ${SLURM_JOB_ID}     # Print final statistics about resource uses before job exits
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
