#!/bin/bash -l
#SBATCH -J kmerGWAS
#SBATCH -e slurm/j%j.run_plink.kmersGWAS.pval_1e-7.err
#SBATCH -o slurm/j%j.run_plink.kmersGWAS.pval_1e-7.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p bigmemm
#SBATCH --time=04-04:00:00
#SBATCH --mem=48G


###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate kmersgwas

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
x=3             # number of fields untils sex (number of _'s +1)
cutoff=34.13    # 0.0000001 p-val chi-squared cutoff (generated from here: http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/pchisq.html)
pval="1e-7"     # must be in #e-# format
# directories
stp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-kmersgwas"
dat_dir="${stp_dir}/data"
out_dir="${stp_dir}/plink_out"
pval_dir="${out_dir}/pval_${pval}"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    working directory : ${dat_dir}"
echo "    output directory  : ${out_dir}" 
echo ""

## ASSOCIATION TEST ##
# navigate to output directory
echo "Checking $out_dir exists...if it doesn't, I guess I'll make it."
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${dat_dir}

# run plink
plink --noweb --bfile sex_kmerGWAS_plink.0 --allow-no-sex --assoc --out ${out_dir}/sex_kmers


## GRAB SIGNIFICANT LOCI ##
# filter large p-values out of association file 
cd ${out_dir}
[ -f sex_kmers_${pval}.assoc.sig ] && rm sex_kmers_${pval}.assoc.sig
awk -v CUTOFF=$cutoff '{ 
if ($8 == "CHISQ") 
    {
    print $0
    }
else if ($8 >= CUTOFF)
    {
    print $0
    }

}' sex_kmers.assoc >> sex_kmers_${pval}.assoc.sig

# sort significant kmers
head -n 1 sex_kmers_${pval}.assoc.sig > sex_kmers_${pval}.assoc.sig.sorted
tail -n +2 sex_kmers_${pval}.assoc.sig | sort -nr -k 8,8 >> sex_kmers_${pval}.assoc.sig.sorted


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
