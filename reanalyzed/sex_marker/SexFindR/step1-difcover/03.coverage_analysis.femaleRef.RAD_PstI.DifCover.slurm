#!/bin/bash -l
#
#SBATCH --array=1-4
#SBATCH -J difCovF
#SBATCH -e slurm/j%j.coverage_analysis.femaleRef.RAD_PstI.DifCover.a%a.err
#SBATCH -o slurm/j%j.coverage_analysis.femaleRef.RAD_PstI.DifCover.a%a.out
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -p bigmemm
#SBATCH --time=02-03:04:05
#SBATCH --mem=16G

# Script to run DifCover on RAD BAM files

###############
###  SETUP  ###
###############
# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate difcover

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
tech="RAD"
ac_line=${SLURM_ARRAY_TASK_ID}
# directories
bam_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step0-mapping/00-alignments_${sex}"
dif_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover"
sex_dir="/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step1-difcover/${SEX}Ref"
sts_dir="${sex_dir}/${tech}_stats"
out_dir="${sex_dir}/${tech}_analysis/a${SLURM_ARRAY_TASK_ID}"
scripts="${dif_dir}/DifCover/dif_cover_scripts"
# files
ac_file="${SEX}Ref_aligned_bam_top2.AC"
###  END CHANGES  ###


#############
###  RUN  ###
#############
# sanity checks
echo ""
echo "Sanity Checks --"
echo "    Sex : ${sex}"
echo "    RAD bam Directory : ${bam_dir}" 
echo "    Output Directory : ${out_dir}"
echo "    DifCover Scripts : ${scripts}"
echo "    Stats Directory : ${sts_dir}"
echo ""

# make sure DifCover has been downloaded
[ -d ${dif_dir} ] || mkdir -p ${dif_dir}
cd ${dif_dir}
[ -d ${scripts} ] || git clone https://github.com/shannonekj/DifCover.git
cd ${scripts}
chmod a+x *
[ -f from_unionbed_to_ratio_per_window_CC0 ] || make

# difcover setup
## assign variables
[ -d ${sts_dir} ] || echo "WARNING: There is no stats directory...backtrack and get some stats."
cd ${sts_dir}
ac=$(sed -n ${ac_line}p ${ac_file} | awk '{print $1}')
fem_bam=$(sed -n ${ac_line}p ${ac_file} | awk '{print $2}')
mal_bam=$(sed -n ${ac_line}p ${ac_file} | awk '{print $3}')
echo "${mal_bam}/${fem_bam} = ${ac}"
## difcover specific variables
FOLDER_PATH=${scripts}
BAM1=${bam_dir}/${fem_bam}
BAM2=${bam_dir}/${mal_bam}
a=10		# minimum coverage for sample1
A=219		# maximum coverage for sample1
b=10		# minimum coverage for sample2
B=240		# maximum coverage for sample2
v=1000		# target number of valid bases in the window
l=500		# minimum size of window to output
AC=${ac}	# Adjustment Coefficient (set 1 if modal coverage is equal) 
p=0.7369656		# enrichment scores threshold (when p=2 will report regions with coverage in sample1 being roughly 4 times larger than coverage in sample2)
bin=1		# for auxiliary stage (5), generates enrichment scores histogram with scores in bins with floating precision 1. For more detailed histogram use 10, 100

# run difcover
[ -d ${out_dir} ] || mkdir -p ${out_dir}
cd ${out_dir}
echo "This directory contains DifCover output from comparing coverage of ${fem_bam} and ${mal_bam} in ${bam_dir}" > README
echo "" >> README
echo "SOFTWARE VERSIONS:" >> README
echo "    -$(samtools --version | head -n 1)" >> README
echo "    -$(samtools --version | head -n 2 | tail -n 1)" >> README
echo "    -$(bedtools --version)" >> README
echo "    -$(R --version | head -n 1)" >> README
echo "" >> README

echo "FILES:" >> README
echo "$(date +%D' '%T) Stage 1:"
call="$FOLDER_PATH/from_bams_to_unionbed.sh $BAM1 $BAM2"
echo "    ${call}"
eval ${call}
cat << EOT >> README
    Renaming.list == file with the original path to sample1 and sample2 sorted bam files
    ref.length.Vk1s_sorted == sorted lengths of reference genome scaffolds
    *.bedcov.Vk1s_sorted == bed coverage files
    *.unionbedcv == coverage of each bed interval for sample1 and sample2 in corresponding columns
EOT
echo "    done."
echo ""

echo "$(date +%D' '%T) Stage 2:"
call="$FOLDER_PATH/from_unionbed_to_ratio_per_window_CC0 sample1_sample2.unionbedcv $a $A $b $B $v $l"
echo "    ${call}"
eval ${call}
echo "    *.ratio_per_w_CC0_a${a}_A${A}_b${b}_B${B}_v${v}_l${l} == ratio of valid bases for sample1 sand sample2 for each window along genome :: calculation ratio=sum_of_cov(sample1)/sum_of_cov(sample2)" >> README
echo "    done."
echo ""

echo "$(date +%D' '%T) Stage 3:"
call="$FOLDER_PATH/from_ratio_per_window__to__DNAcopy_output.sh "sample1_sample2.ratio_per_w_CC0_a"$a"_A"$A"_b"$b"_B"$B"_v"$v"_l"$l $AC"
echo "    ${call}"
eval ${call}
cat << EOT >> README
    *.log2adj_${AC} == AC*log2(ratio) of above file for each window
    *.DNAcopyout == [chr | loc.start | loc.end | num.mark | seg.mean]; num.mark=# markers in segment; seg.mean=mean copy number ratio of each segment
    *.pdf == plots of whole and plateau copy number ratios between sample1 and sample2; plot one is ordered by locations throughout genome; plot two is ordered by the value of the segment means
EOT
echo "    done."
echo ""

echo "$(date +%D' '%T) Stage 4:"
call="$FOLDER_PATH/get_DNAcopyout_with_length_of_intervals.sh "sample1_sample2.ratio_per_w_CC0_a"$a"_A"$A"_b"$b"_B"$B"_v"$v"_l"$l".log2adj_"$AC".DNAcopyout" ref.length.Vk1s_sorted"
echo "    ${call}"
eval ${call}
echo "    *.DNAcopyout.len == [Contig | contig_length | begin_of_interval|  length_of_interval | ratio]; length of each DNAcopy intervals with corresponding ratios" >> README
call="$FOLDER_PATH/generate_DNAcopyout_len_histogram.sh "sample1_sample2.ratio_per_w_CC0_a"$a"_A"$A"_b"$b"_B"$B"_v"$v"_l"$l".log2adj_"$AC".DNAcopyout.len" $bin"
echo "    ${call}"
eval ${call}
echo "    *.DNAcopyout.len.hist_b${bin} == length of each DNA copy interval with corresponding whole contig (DifCover explanation = computed length of each DNAcopy interval and output length of the whole contig to which that interval belongs)" >> README
echo "    done."
echo ""

echo "$(date +%D' '%T) Stage 5:"
call="$FOLDER_PATH/from_DNAcopyout_to_p_fragments.sh "sample1_sample2.ratio_per_w_CC0_a"$a"_A"$A"_b"$b"_B"$B"_v"$v"_l"$l".log2adj_"$AC".DNAcopyout" $p"
echo "    ${call}"
eval ${call}
cat << EOT >> README
    *.DNAcopyout.up${p} == sample1_coverage > sample2_coverage; DNAcopy intervals with av(adj_coef*log2ratio)> p, where p=enrichment scores threshold (when p=2 will report regions with coverage in sample1 being roughly 4 times larger than coverage in sample2)
    *.DNAcopyout.down${p} == sample2_coverage > sample1_coverage; DNAcopy intervals with av(adj_coef*log2ratio)> p, where p=enrichment scores threshold (when p=2 will report regions with coverage in sample1 being roughly 4 times larger than coverage in sample2)
EOT
echo "    done."
echo ""


# print job information
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
end=`date +%s`
runtime=$((end-start))
echo "Runtime : "${runtime}
