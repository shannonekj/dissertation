#!/bin/bash -l
#SBATCH -J cntKmers
#SBATCH -e slurm/j%j.count_kmers.kmersGWAS.err
#SBATCH -o slurm/j%j.count_kmers.kmersGWAS.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=2G


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
n=5 # number of _'s before R1 in fastq file names
# directories
dat_dir="/group/millermrgrp2/shannon/projects/DS_sex-marker/00-RAD-data"
stp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-kmersgwas"
out_dir="${stp_dir}/data"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    fastq directory : ${dat_dir}" 
echo "    output directory : ${out_dir}"
echo ""

# navigate to output directory
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}

# setup directory structure
# make directories for fastq's
for i in $(ls --color=none ${dat_dir}/*_R1.fastq)
do 
    baseName=$(echo ${i} | rev | cut -d/ -f1 | rev | cut -d_ -f1-${n})
    [ -d ${baseName} ] || mkdir -p ${baseName}
done


## create lists of input fastq's to count kmers
ls --color=none -d */ > dirlist.txt
sed 's/\///g' dirlist.txt > clean_dirlist.txt
for dir in $(cat clean_dirlist.txt)
do
    cd ${dir}
    echo ${dir}
    ls --color=none ${dat_dir}/${dir}*.fastq > input_files.txt
    cd ..
done

# count kmers 
for dir in $(cat clean_dirlist.txt)
do 
    cd ${dir}
    cat << EOT >> step1_kmerGWAS.sh
#!/bin/bash
#
#SBATCH -J canon_KWAS 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 16
#SBATCH -p bigmemm
#SBATCH -t 0-11:59:00
#SBATCH --mem=64G                 # Memory pool for all cores 
#SBATCH -o CanonKWAS.%A.out
#SBATCH -e CanonKWAS.%A.err

echo \$\(date +%D' '%T\) "Counting canonical kmers occuring more than n=2"
${stp_dir}/kmersGWAS/external_programs/kmc_v3 -t16 -k31 -ci2 @input_files.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
echo "    done."
echo ""
echo \$\(date +%D' '%T\) "Counting all kmers w/o canonical transformation"
${stp_dir}/kmersGWAS/external_programs/kmc_v3 -t16 -k31 -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2
echo "    done."
echo ""
echo \$\(date +%D' '%T\) "Adding strand information"
${stp_dir}/kmersGWAS/bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand
echo "    done."
EOT
    sbatch step1_kmerGWAS.sh
    cd ../
done

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
