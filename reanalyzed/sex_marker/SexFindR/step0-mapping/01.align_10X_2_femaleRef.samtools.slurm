#!/bin/bash -l
#SBATCH -J aln2fm
#SBATCH -e slurm/j%j.align_10X_2_femaleRef.samtools.err
#SBATCH -o slurm/j%j.align_10X_2_femaleRef.samtools.out
#SBATCH --ntasks=16
#SBATCH -p bigmemh
#SBATCH --time=5-04:03:02
#SBATCH --mem=24G

# Script to align 10X data to female reference genome using bwa & samtools
# SOFTWARE VERSIONS:
##      -samtools=1.10
##      -bwa 0.7.17-r1188

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
samtools --version

### MAKE CHANGES ###
sex="F"
SEX="female"
# directories
tenx_dir="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_${sex}"
out_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-mapping/10X-alignments_${sex}"
# files
if [ $sex == "M" ]
then
    reference="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/04-linkage_map/out_chromonomer/CHRR_integrated.fa"
    r1="Male2_S63_L004_R1_001.fastq.gz"
    r2="Male2_S63_L004_R2_001.fastq.gz"
elif [ $sex == "F" ]
then
    reference="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/03-assemblies/sandbox_ipa_F/Hyp_tra_F_combined_chromonomer/output/CHRR_integrated.fa"
    r1="Smelt2F_S432_L004_R1_001.fastq.gz"
    r2="Smelt2F_S432_L004_R2_001.fastq.gz"
else
    echo "sex input is not female or male. exiting process."
    exit
fi
###  END CHANGES  ###


#############
###  RUN  ###
#############

# make output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# index reference genome
echo "$(date +%D' '%T) indexing reference genome"
## symlink reference genome
ref_short=$(echo $reference | rev | cut -f1 -d/ | rev)
ref_pre=$(echo $ref_short | cut -f1 -d.)_${sex}.fa
## make symbolic link to reference genome
[ -d ${out_dir}/ref_genome ] || mkdir -p ${out_dir}/ref_genome
cd ${out_dir}/ref_genome
[ -h ${ref_pre} ] || ln -s ${reference} ${ref_pre}
## index reference genome
[ -f ${ref_pre}.sa ] || bwa index ${ref_pre}; sleep 1m
echo "    done."

# align files
cd ${out_dir}
echo "$(date +%D' '%T) aligning ${SEX} 10X files"  
bwa mem -t ${SLURM_NTASKS} ${out_dir}/ref_genome/${ref_pre} ${tenx_dir}/${r1} ${tenx_dir}/${r2} | samtools view -Sb - | samtools sort -n - | samtools fixmate -m - tenx_${sex}.sort-n.fixmate-m.bam
echo "    done."
echo ""

echo "$(date +%D' '%T)  Sorting and removing duplicates for ${SEX} 10X data"
samtools sort tenx_${sex}.sort-n.fixmate-m.bam | samtools markdup -r - tenx_${sex}.sort-n.fixmate-m.sort.markdup-r.bam
echo "    done."
echo ""

echo $(date +%D' '%T)  Indexing bam file
samtools index tenx_${sex}.sort-n.fixmate-m.sort.markdup-r.bam
echo "    done."
echo ""


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
