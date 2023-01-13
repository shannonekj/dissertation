#!/bin/bash
#
#SBATCH -J indVCF
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-02:59:59
#SBATCH --mem=4Gb
#SBATCH -o slurm/j%j.get_SNPdensity.out
#SBATCH -e slurm/j%j.get_SNPdensity.err

# $1 is name of file containing sample list
# $2 is output vcf name
# $3 is input vcf name
# $4 is sex

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate bcftools

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

echo SCRIPT FOR : $1

echo "Step 1 : creating individual vcf file"
bcftools view -a -s $1 -o individual_vcfs/$4_$2 $3

echo "Step 2 : compressing"
bgzip -c individual_vcfs/$4_$2 > individual_vcfs/$4_${2}.gz

echo "Step 3 : indexing"
bcftools index individual_vcfs/$4_${2}.gz

echo "Step 4 : getting SNP density"
vcftools --vcf individual_vcfs/$4_$2 --SNPdensity 10000 --out $4_$2
