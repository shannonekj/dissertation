#!/bin/bash
#SBATCH --array=1-2
#SBATCH -J srmsh_sg
#SBATCH -e slurm/06.sourmash_signatures.j%j.err
#SBATCH -o slurm/06.sourmash_signatures.j%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=10-24:24:24
#SBATCH -p bigmemh


# commands for kmer analysis
## setup environment
### initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

### activate conda env
conda activate sourmash

### fail on weird errors
set -o nounset
set -o errexit
set -x

### directories
wrk_dir="/group/millermrgrp2/shannon/projects/DS_sex-marker/04-sourmash_indiv"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/kmer_analysis"


### array things
files=(Smelt Male)
sexes=(female male)
x=${SLURM_ARRAY_TASK_ID}
echo "x is ${x}"
prefix=${files[$(( $x - 1 ))]}
sex=${sexes[$(( $x - 1 ))]}

echo "file prefix is $prefix"
echo "sex is $sex"

### make output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

## run sourmash sig
echo $(date +%D' '%T) "Running sourmash sig merge..."
### sourmash sig merge - merge R1 & R2 sigs for each sex
sourmash sig merge ${wrk_dir}/${prefix}*.sample.sig -o ${sex}.sig -k 21
echo $(date +%D' '%T) "    Success!"
    # output = {sex}.sig
echo $(date +%D' '%T) "Running sourmash sig filter..."
### sourmash sig filter - filter on k-mer abundance of 5000
sourmash sig filter ${sex}.sig -m 5 -o ${sex}-1k-m5.sig
echo $(date +%D' '%T) "    Success!"
    # output = {sex}-1k-m5.sig

# now do do things in python
