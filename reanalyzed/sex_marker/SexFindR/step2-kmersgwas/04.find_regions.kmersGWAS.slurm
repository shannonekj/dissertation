#!/bin/bash -l
#SBATCH -J kmerGWAS
#SBATCH -e slurm/j%j.find_regions.kmersGWAS.pval_1e-5.err
#SBATCH -o slurm/j%j.find_regions.kmersGWAS.pval_1e-5.out
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
pval="1e-5"       # must be in #e-# format
cutoff=25.0135    # 0.00001 p-val chi-squared cutoff (generated from here: http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/pchisq.html)
prefix="sex_kmers_${pval}"
# directories
stp_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-kmersgwas"
plk_dir="${stp_dir}/plink_out"
pval_dir="${plk_dir}/pval_${pval}"
aby_dir="${stp_dir}/ABYSS_out/pval_${pval}"
bst_dir="${stp_dir}/blast_out/pval_${pval}"
# files
fem_ref="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/03-assemblies/sandbox_ipa_F/Hyp_tra_F_combined_chromonomer/output/CHRR_integrated.fa"
mal_ref="/group/millermrgrp2/shannon/projects/assembly_genome_Hypomesus-transpacificus/04-linkage_map/out_chromonomer/CHRR_integrated.fa"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    plink directory : ${plk_dir}"
echo "    p-value directory : ${pval_dir}"
echo "    ABYSS directory : ${aby_dir}"
echo "    blast directory : ${bst_dir}"
echo "    female ref genome : ${fem_ref}"
echo "    male ref genome : ${mal_ref}"
echo "    chi-square cutoff : ${cutoff}"
echo ""

# navigate to output directory
echo "Checking all the directories exist...if they don't, I guess I'll make them."
[ -d ${pval_dir} ] || mkdir -p ${pval_dir}
[ -d ${aby_dir} ] || mkdir -p ${aby_dir}
[ -d ${bst_dir} ] || mkdir -p ${bst_dir}
echo ""

# grab significant loci
## filter
cd ${plk_dir}
[ -f ${pval_dir}/sex_kmers_${pval}.assoc.sig ] && rm ${pval_dir}/sex_kmers_${pval}.assoc.sig
awk -v CUTOFF=$cutoff '{
if ($8 == "CHISQ")
    {
    print $0
    }
else if ($8 >= CUTOFF)
    {
    print $0
    }
}' sex_kmers.assoc >> ${pval_dir}/sex_kmers_${pval}.assoc.sig

## sort 
cd ${pval_dir}
head -n 1 sex_kmers_${pval}.assoc.sig > sex_kmers_${pval}.assoc.sig.sorted
tail -n +2 sex_kmers_${pval}.assoc.sig | sort -nr -k 8,8 >> sex_kmers_${pval}.assoc.sig.sorted


# ABYSS
cd ${aby_dir}
## format input file
echo $(date +%D' '%T) "Formatting files for ABYSS."
[ -h ${prefix}.assoc.sig ] || ln -s ${pval_dir}/${prefix}.assoc.sig ${prefix}.assoc.sig
tail -n +2 ${prefix}.assoc.sig | awk 'BEGIN { OFS="" } { print "\>",$9, "\n",$2 }' > ${prefix}_abyss_input.txt
echo "    done."
echo ""

## run ABYSS
echo $(date +%D' '%T) "Running ABYSS"
ABYSS -k25 -c0 -e0 ${prefix}_abyss_input.txt -o ${prefix}_abyss_output.txt
echo "    done."
echo ""


# BLAST
cd ${bst_dir}
## female reference
echo $(date +%D' '%T) "Locating k-mers in female reference genome"
### symlink
echo $(date +%D' '%T) "Symlinking reference"
sex="F"
SEX="female"
reference=${fem_ref}
ref_short=$(echo $reference | rev | cut -f1 -d/ | rev)
ref_pre=$(echo $ref_short | cut -f1 -d.)_${sex}.fa
mkdir -p ${SEX}Ref
cd ${SEX}Ref
[ -h ${ref_pre} ] || ln -s ${reference} ${ref_pre}
echo "    done."
echo ""
### make blastdb
echo $(date +%D' '%T) "Making db"
makeblastdb -in ${ref_pre} -dbtype nucl
echo "   done."
echo ""
### run BLAST
echo $(date +%D' '%T) "BLASTing"
blastn -query ${aby_dir}/${prefix}_abyss_output.txt -db ${ref_pre} -outfmt 6 >> ${prefix}_blast_results.${SEX}Ref.outfmt6
echo "    done."
### sort results
echo $(date +%D' '%T) "Sorting BLAST results"
sort -k2,2 ${prefix}_blast_results.${SEX}Ref.outfmt6 >> ${prefix}_blast_results.${SEX}Ref.outfmt6.sorted
echo "    done."
echo ""

## male reference
cd ${bst_dir}
echo $(date +%D' '%T) "Locating k-mers in male reference genome"
### symlink
echo $(date +%D' '%T) "Symlinking reference"
sex="M"
SEX="male"
reference=${mal_ref}
ref_short=$(echo $reference | rev | cut -f1 -d/ | rev)
ref_pre=$(echo $ref_short | cut -f1 -d.)_${sex}.fa
mkdir -p ${SEX}Ref
cd ${SEX}Ref
[ -h ${ref_pre} ] || ln -s ${reference} ${ref_pre}
echo "    done."
echo ""
### make blastdb
echo $(date +%D' '%T) "Making db for running local alignment"
makeblastdb -in ${ref_pre} -dbtype nucl
echo "   done."
echo ""
### run BLAST
echo $(date +%D' '%T) "BLASTing"
blastn -query ${aby_dir}/${prefix}_abyss_output.txt -db ${ref_pre} -outfmt 6 >> ${prefix}_blast_results.${SEX}Ref.outfmt6
echo "    done."
echo  ""
### sort results
echo $(date +%D' '%T) "Sorting BLAST results"
sort -k2,2 ${prefix}_blast_results.${SEX}Ref.outfmt6 >> ${prefix}_blast_results.${SEX}Ref.outfmt6.sorted
echo "    done."
echo ""

# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
