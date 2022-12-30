#!/bin/R
# Taken from https://github.com/drtamermansour/Equine_parentage_testing/blob/main/scripts/rainfall.R


# SETUP
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16") # note you need to be running R 4.3 (r-devel) version
#BiocManager::install("karyoploteR")
#BiocManager::install("Hmisc")
library(karyoploteR)
library(tidyverse)


# INPUT : BED file
### chr   start   end
### lg01  1       154687332
### lg02  1       74835431
### ...


# READ IN DATA
setwd("/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis")
## genome
male_genome <- read.table(file="CHRR_integrated.male.fa.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
male_genome <- setNames(male_genome, c("chr", "start", "end"))
## putative Y
put_Y <- read.table(file="putative_y.maleRef.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
put_Y <- setNames(put_Y, c("chr", "start_pos", "end_pos", "read_name", "score", "strand"))
pY_regs <-toGRanges(put_Y)
### extract unique chromosomes which putative Y seq mapped to
#Y_chrs <- unique(put_Y$chr)
#regions <- c()
#for (x in 1:length(Y_chrs)) {
#  chromo <- Y_chrs[x]
#  if (x == 1) {
#    regions <- paste('"', chromo, '"', sep="")
#  } else {
#    regions <- paste(regions, ', "', chromo, '"', sep="")
#  }
#}


# PLOT STUFF
## All Chromosomes
### filtered  to a reasonable number of scaffolds (n=value(mean))
mg_filt <- male_genome %>% filter(end >= mean(end, na.rm=TRUE))
mg_filt <- toGRanges(data.frame(chr=mg_filt$chr, start=mg_filt$start, end=mg_filt$end))
mg_karyotype <- plotKaryotype(genome=mg_filt, plot.type=6)
## Just Chromosomes where putative Y mapped
#### MANUAL ENTRY <- I suck
kp <- plotKaryotype(genome=male_genome, chromosomes=c("scaffold_9", "lg09", "lg17", "lg21", "scaffold_190", "lg23"))
kpPlotRegions(kp, data=pY_regs)
kpAddBaseNumbers(kp, tick.dist=1000000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=100000, minor.tick.len=5, minor.tick.col="gray")
## Zoomed
# I could have a facet plot that has 6 regions, then have a for loop that iterates to plot each manually based on # in vector... 
### scaffold_9
#### extract scaffold from put_y to get range
scf9 <- "scaffold_9"
scf9_region <- put_Y %>% filter(chr == scf9, na.rm=TRUE)
min(scf9_region$end_pos)
max(scf9_region$start_pos)
#### get bams
scf9_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.scaffold_9.merged.subset.bam"
scf9_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.scaffold_9.merged.subset.bam"
#### define zoom area
#zoom.region <- toGRanges(data.frame(genome=male_genome, chromosomes="scaffold_9", beg_reg, end_reg))
scf9_kp <- plotKaryotype(genome=male_genome, chromosomes="scaffold_9", zoom="scaffold_9:6110000-6420000")
kpPlotRegions(scf9_kp, data=pY_regs)
kpAddBaseNumbers(scf9_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(scf9_kp, data=scf9_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(scf9_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(scf9_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(scf9_kp, data=scf9_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(scf9_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(scf9_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)

### "lg09"
lg09 <- "lg09"
lg09_region <- put_Y %>% filter(chr == lg09, na.rm=TRUE)
min(lg09_region$start_pos)
max(lg09_region$end_pos)
lg09_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.lg09.merged.subset.bam"
lg09_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.lg09.merged.subset.bam"
lg09_kp <- plotKaryotype(genome=male_genome, chromosomes="lg09", zoom="lg09:21000000-21100000")
kpPlotRegions(lg09_kp, data=pY_regs)
kpAddBaseNumbers(lg09_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(lg09_kp, data=lg09_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(lg09_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(lg09_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(lg09_kp, data=lg09_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(lg09_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(lg09_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)

### "lg17"
lg17_region <- put_Y %>% filter(chr == "lg17", na.rm=TRUE)
min(lg17_region$start_pos)
max(lg17_region$end_pos)
lg17_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.lg17.merged.subset.bam"
lg17_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.lg17.merged.subset.bam"
lg17_kp <- plotKaryotype(genome=male_genome, chromosomes="lg17", zoom="lg17:8100000-8200000")
kpPlotRegions(lg17_kp, data=pY_regs)
kpAddBaseNumbers(lg17_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(lg17_kp, data=lg17_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(lg17_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(lg17_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(lg17_kp, data=lg17_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(lg17_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(lg17_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)

### "lg21"
lg21_region <- put_Y %>% filter(chr == "lg21", na.rm=TRUE)
min(lg21_region$start_pos)
max(lg21_region$end_pos)
lg21_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.lg21.merged.subset.bam"
lg21_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.lg21.merged.subset.bam"
lg21_kp <- plotKaryotype(genome=male_genome, chromosomes="lg21", zoom="lg21:9815000-9860000")
kpPlotRegions(lg21_kp, data=pY_regs)
kpAddBaseNumbers(lg21_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(lg21_kp, data=lg21_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(lg21_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(lg21_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(lg21_kp, data=lg21_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(lg21_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(lg21_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)

### "lg23"
lg23_region <- put_Y %>% filter(chr == "lg23", na.rm=TRUE)
min(lg23_region$start_pos)
max(lg23_region$end_pos)
lg23_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.lg23.merged.subset.bam"
lg23_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.lg23.merged.subset.bam"
lg23_kp <- plotKaryotype(genome=male_genome, chromosomes="lg23", zoom="lg23:10550000-10620000")
kpPlotRegions(lg23_kp, data=pY_regs)
kpAddBaseNumbers(lg23_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(lg23_kp, data=lg23_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(lg23_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(lg23_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(lg23_kp, data=lg23_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(lg23_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(lg23_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)

### "scaffold_190"
scf190_region <- put_Y %>% filter(chr == "scaffold_190", na.rm=TRUE)
min(scf190_region$start_pos)
max(scf190_region$end_pos)
scf190_female <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/female.scaffold_190.merged.subset.bam"
scf190_male <- "/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/kmer_analysis/male.scaffold_190.merged.subset.bam"
scf190_kp <- plotKaryotype(genome=male_genome, chromosomes="scaffold_190", zoom="scaffold_190:80000-100000")
kpPlotRegions(scf190_kp, data=pY_regs)
kpAddBaseNumbers(scf190_kp, tick.dist=10000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=1000, minor.tick.len=5, minor.tick.col="gray")
kpPlotBAMCoverage(scf190_kp, data=scf190_male, col="#1E90FF", border=NA, r0=0.5, r1=0, ymax=600)
kpAxis(scf190_kp, r0=0.5, r1=0, ymax=600)
kpAddLabels(scf190_kp, "Male", r0=0, r1=0.5, label.margin = 0.05)
kpPlotBAMCoverage(scf190_kp, data=scf190_female, col="#9666FF", border=NA, r0=0.5, r1=1, ymax=600)
kpAxis(scf190_kp, r0=0.5, r1=1, ymax=600)
kpAddLabels(scf190_kp, "Female", r0=0.5, r1=1, label.margin = 0.05)
