#!/bin/bash -l
#SBATCH -J vcFem
#SBATCH -e slurm/j%j.variant_calling.femaleRef.platypus.err
#SBATCH -o slurm/j%j.variant_calling.femaleRef.platypus.out
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -p bigmemh
#SBATCH --time=24:00:00
#SBATCH --mem=16G

# Script to call variants from RAD reads aligned w/samtools using female reference genome

###############
###  SETUP  ###
###############

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate platypus

# fail on weird errors
set -o nounset
set -o errexit
set -x
set -e

hostname
start=`date +%s`
echo "My SLURM_JOB_ID: $SLURM_JOB_ID"

### MAKE CHANGES ###
sex="F"
SEX="female"
bam_suffix=".sort-n.fixmate-m.sort.markdup-r.bam"
bamlist="${SEX}Ref.bamlist"
# directories
bam_dir="/home/sejoslin/git_repos/SexFindR/analysis_2022/step0-mapping/00-alignments_${sex}"
out_dir="/home/sejoslin/git_repos/SexFindR/analysis_2022/step0-variantcalling/${SEX}Ref"
plat_dir="/home/sejoslin/miniconda3/envs/platypus/share/platypus-variant-0.8.1.2-0"
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
echo "Sex : ${sex}"
echo "Reference Genome : ${reference}"

echo $(date +%D' '%T) "Setting up directories and symlinks"
# make output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}
# symlink reference genome
ref_short=$(echo $reference | rev | cut -f1 -d/ | rev)
ref_pre=$(echo $ref_short | cut -f1 -d.)_${sex}
[ -h ${ref_short} ] || ln -s ${reference} ${ref_short}
echo "    done."
echo ""

# make bamlist of all RAD seq files
echo $(date +%D' '%T) "Making bamlist"
ls ${bam_dir}/*${bam_suffix} > ${bamlist}
echo "    done."
echo ""

# index referenc
echo $(date +%D' '%T) "Indexing reference"
conda activate samtools
samtools faidx ${ref_short}
echo "    done."
echo ""

# call variants
echo $(date +%D' '%T) "Calling Variants"
conda activate platypus
python ${plat_dir}/Platypus.py callVariants --verbosity=3 --bufferSize=50000 --nCPU=${SLURM_NTASKS} --bamFiles=${bamlist} --refFile=${ref_short} --output=${SEX}Ref_RAD-Pst1_24F_24M.vcf
echo "    done."

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
