#!/bin/R

# USAGE: Rscript plot_difcover.R ${DNAcopyout} ${scaffold_length_file} ${out_dir}

options(repos=c(CRAN="https://ftp.osuosl.org/pub/cran/"))

# define variables
args<-commandArgs(TRUE)
DNA_out <- args[1]
scaf_file <- args[2]
out_dir <- args[3]
P <- args[4]
a <- args[5]
p <- as.numeric(P)

# sanity check
print(paste("DNAcopy output file:", DNA_out, sep=" "))
print(paste("Scaffold Lengths file:", scaf_file, sep=" "))
print(paste("Output Directory:", out_dir, sep=" "))
print(paste("enrichment scores threshold:", p, sep=" "))
str(p)
print(paste("array:", a, sep=" "))

# require packages
list.of.packages <- c("tidyverse", "patchwork", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(tidyverse)
library(patchwork)
library(ggpubr)

setwd(out_dir)
getwd()
together <- read_tsv(file=DNA_out, col_names=F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2.MaleCoverage_over_FemaleCoverage"=X5) %>% mutate("bases spanned" = stop-base)
head(together)
scaffold_lengths <- read_tsv(scaf_file, col_names=c("scaf", "length"))
head(scaffold_lengths)
proportion <- full_join(together, scaffold_lengths) %>% mutate(Chromosome = ifelse(grepl("scaffold_", scaf), "unplaced", "Autosome")) %>% mutate(proportion = `bases spanned`/length)
head(proportion)

s <- 35
s

x <- ggscatter(proportion,
               x = "log2.MaleCoverage_over_FemaleCoverage",
               y = "proportion",
               color = "Chromosome",
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"),
               ylab = "Proportion of chromosome (in region)",
               size = 2)  + geom_vline(xintercept =0,linetype="dotted")  +
                font("legend.title", size = s) +
                font("legend.text", size = s) +
                font("xlab", size = s) +
                font("ylab", size = s) +
                font("xy.text", size = s)

filtered_proportion <- proportion %>% filter(log2.MaleCoverage_over_FemaleCoverage >= p | log2.MaleCoverage_over_FemaleCoverage <= -p) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)
head(filtered_proportion)

y <- ggdotchart(filtered_proportion,
                x = "scaf",
                y = "total chromosome proportion with significantly different coverage",
                color = "scaf",
                #palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"),
                xlab = "Chromosome",
                ylab = "Proportion of chromosome (significant)", 
                #sorting = "descending",
                add = "segments",
                add.params = list(color = "lightgray", size = 1),
                group = "Chromosome",
                dot.size = 4 )  +
                font("legend.title", size = s) +
                font("legend.text", size = s) +
                font("xlab", size = s) +
                font("ylab", size = s) +
                #theme(axis.text.x=element_blank()) +
                font("xy.text", size = s) 


pdf(paste("DNAcopyout.a", a, ".pdf", sep=""), width = 23, height = 11.13349)
z <- x+y
print(z)
dev.off()
dev.off()
#ggsave(paste(a, DNA_out, "pdf", sep="."), width = 23, height = 11.13349)

