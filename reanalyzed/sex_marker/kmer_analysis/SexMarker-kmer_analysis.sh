#!/bin/bash
#SBATCH --array=1-4
#SBATCH -J srmsh_sx
#SBATCH -e slurm/06.sourmash_signatures.j%j.err
#SBATCH -o slurm/06.sourmash_signatures.j%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=10-24:24:24
#SBATCH -p high


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
tenx_m_dir="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_M"
tenx_m_r1="Male2_S63_L004_R1_001.fastq.gz"
tenx_m_r2="Male2_S63_L004_R2_001.fastq.gz"
tenx_f_dir="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_F"
tenx_f_r1="Smelt2F_S432_L004_R1_001.fastq.gz"
tenx_f_r2="Smelt2F_S432_L004_R2_001.fastq.gz"

out_dir="/group/millermrgrp2/shannon/projects/DS_sex-marker/04-sourmash_indiv"

### array things
files=(/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_M/Male2_S63_L004_R1_001.fastq.gz /group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_M/Male2_S63_L004_R2_001.fastq.gz /group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_F/Smelt2F_S432_L004_R1_001.fastq.gz /group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/00-raw_data/data-10X_F/Smelt2F_S432_L004_R2_001.fastq.gz)
x=${SLURM_ARRAY_TASK_ID}
echo "x is ${x}"
f=${files[$(( $x - 1 ))]}
echo "Global path to fastq is $f"
file=$(echo $f | rev | cut -d/ -f1 | rev)
echo "file name is $file"

### make output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

### sym link tenx files
ln -s ${f} ${file}


## run sourmash
### sourmash compute - create MinHash sketches at k-mer size 
echo $(date +%D' '%T) "starting sourmash compute"
sourmash compute -k 21,31,51 --scaled 1000 --track-abundance ${file} -o ${file}.sample.sig
echo $(date +%D' '%T) "    sourmash run complete!"
    # output = {R1_file}.sample.sig & {R2_file}.sample.sig for each sex
### sourmash sig merge - merge R1 & R2 sigs for each sex
sourmash sig merge Smelt2* -o female.sig -k 21
sourmash sig merge Male2* -o male.sig -k 21
    # output = {sex}.sig
### sourmash sig filter - filter on k-mer abundance of 5000
sourmash sig filter female.sig -m 5 -o female-1k-m5.sig
sourmash sig filter male.sig -m 5 -o male-1k-m5.sig

# now do do things in python
