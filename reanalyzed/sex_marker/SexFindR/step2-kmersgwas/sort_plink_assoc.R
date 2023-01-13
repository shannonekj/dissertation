#!/bin/R

# USAGE: Rscript sort_plink_assoc.R ${association_file} ${out_dir}

args <- commandArgs(TRUE)
infile <- args[1]
outdir <- args[2]
outfile <- paste(infile, "sig", sep=".")

#sanity checks
print(paste("Association file :", infile, sep=" "))
print(paste("Output Directory:", outdir, sep=" "))
print(paste("Output file :", outfile, sep=" "))

# setup env
library(tidyverse)
setwd(outdir)
getwd()
assoc <- read.table(infile, header=T)


#save
write.table(sig, file=paste(infile, "sig", sep=".") 
