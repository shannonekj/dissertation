#!/bin/bash -l
#SBATCH -J combineK
#SBATCH -e slurm/j%j.combine_kmers.kmersGWAS.err
#SBATCH -o slurm/j%j.combine_kmers.kmersGWAS.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p bigmemm
#SBATCH --time=04:00:00
#SBATCH --mem=12G


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
x=3     # number of fields untils sex (number of _'s +1)
# directories
stp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-kmersgwas"
out_dir="${stp_dir}/data"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    working directory : ${out_dir}"
echo ""

# navigate to output directory
echo "Checking $out_dir exists..."
[ -d ${out_dir} ] || exit
echo "    it does."
cd ${out_dir}

# make list of kmers_with_strand locations
## File FMT: [ /global/path/to/individual_id/kmers_with_strand      individual_id
echo $(date +%D' '%T) "Making a list of kmer with strand files"
[ -f final_kmers_clean_list_path.txt ] && rm final_kmers_clean_list_path.txt
for i in $(cat clean_dirlist.txt)
do
    echo -e "${PWD}/${i}/kmers_with_strand\t${i}" >> final_kmers_clean_list_path.txt
done
echo "   done."
echo ""

# combine kmers
echo $(date +%D' '%T) "Combining kmers"
${stp_dir}/kmersGWAS/bin/list_kmers_found_in_multiple_samples -l final_kmers_clean_list_path.txt -k 31 --mac 5 -p 0.2 -o kmers_to_use
echo "    done."
echo ""

# create kmer tables
echo $(date +%D' '%T) "Creating kmer tables"
${stp_dir}/kmersGWAS/bin/build_kmers_table -l final_kmers_clean_list_path.txt -k 31 -a kmers_to_use -o kmers_table
echo "    done."
echo ""

# make plink files
## phenotype.pheno
echo $(date +%D' '%T) "Making PHENO file"
echo -e "accession_id\tphenotype_value" > phenotype.pheno
for i in $(cat clean_dirlist.txt)
do
    sex=$(echo ${i} | cut -d_ -f${x})
    if [ $sex == "F" ]
    then
        val="1"
    elif [ $sex == "M" ]
    then
        val="2"
    else
        echo "sex input is not female or male. exiting process."
        exit
    fi
    echo -e "${i}\t${val}" >> phenotype.pheno
done
echo "    done."
echo ""

## bed files
echo $(date +%D' '%T) "Making BED files"
${stp_dir}/kmersGWAS/bin/kmers_table_to_bed -t kmers_table -k 31 -p phenotype.pheno --maf 0.05 --mac 5 -b 1000000000 -o sex_kmerGWAS_plink
echo "    done."


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
