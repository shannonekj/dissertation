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

pY_regs <- toGRanges(data.frame(chr=put_Y$chr, start=put_Y$start_pos, end=put_Y$end_pos)
toGRanges(data.frame(chr=mg_filt$chr, start=mg_filt$start, end=mg_filt$end))
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
### plot
#### MANUAL ENTRY <- I suck
kp <- plotKaryotype(genome=male_genome, chromosomes=c("scaffold_9", "lg09", "lg17", "lg21", "scaffold_190", "lg23"))


Y_regs <- toGRanges(c(regions))




## filter to a reasonable number of scaffolds
mg_filt <- male_genome %>% filter(end >= mean(end, na.rm=TRUE))
mg_filt <- toGRanges(data.frame(chr=mg_filt$chr, start=mg_filt$start, end=mg_filt$end))
mg_karyotype <- plotKaryotype(genome=mg_filt, plot.type=6)

## testing things on Chr 8
regs <- toGRanges(c("lg08:1000000-2000000", "lg08:300000-500000", "lg08:7000000-8500000"))
colors <- c("red", "#889F34", lighter(rainbow(n = 18)[12], 50))
mg_eight <- plotKaryotype(genome=mg_filt, chromosomes="lg08")
kpPlotRegions(mg_eight, data=regs, r0=0, r1=0.45, col = lighter(colors), border=colors, lwd=3)
kpAddBaseNumbers(mg_eight, tick.dist=1000000, tick.len=10, tick.col="red", cex=1, 
                 minor.tick.dist=100000, minor.tick.len=5, minor.tick.col="gray")

## Load input VCF and tranform into GRange object
cand.snps=read.table(file="Equ_Parentv2cor.check.vcf", header=FALSE, sep="\t", stringsAsFactors=FALSE)
cand.snps <- setNames(cand.snps, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
cand.snps[cand.snps == "eMSYv3"] <- "Y"
library(regioneR)
cs.gr <- toGRanges(cand.snps[,c("CHROM", "POS", "POS", "ID", "REF", "ALT")])

## generate custom genome (as a GRange object as well)
equ_genome=grep('contig',readLines("Equ_Parentv2cor.check.vcf"), value = TRUE)
equ_genome=gsub(",|>","=",equ_genome)
equ_genome<-read.table(textConnection(equ_genome),comment.char = "",sep = "=")
equ_genome[,4]=1
equ_genome[equ_genome == "eMSYv3"] <- "Y"
equ_genome <- toGRanges(data.frame(chr=equ_genome$V3,start=equ_genome$V4,end=equ_genome$V5))

## Starting the ploting
library(karyoploteR)

## Create a simple Rainfall plot
kp <- plotKaryotype(plot.type=4,genome=equ_genome, chromosomes=c(1:31,"X","Y"))
kpPlotRainfall(kp, data = cs.gr)

## Create a simple Rainfall plot but control the colors and add labels
#variant.colors <- getVariantsColors(cs.gr$REF, cs.gr$ALT)
#kp <- plotKaryotype(plot.type=4,genome=equ_genome, chromosomes=c(1:31,"X","Y"))
#kpPlotRainfall(kp, data = cs.gr, col=variant.colors)
#kpAddLabels(kp, labels = c("Distance between markers (log10)"), srt=90, pos=1, label.margin = 0.04)

## Create a Rainfall plot with custom plotting parameters
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20
pp$leftmargin <- 0.08
#kp <- plotKaryotype(plot.type=4, genome=equ_genome, chromosomes=c(1:31,"X","Y"),
#                    ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
kp <- plotKaryotype(plot.type=4, genome=equ_genome, chromosomes=c(1:31,"X","Y"),
                    labels.plotter = NULL, plot.params = pp)
#kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=0.7)

#kpPlotRainfall(kp, data = cs.gr, col=variant.colors)
#kpPlotRainfall(kp, data = cs.gr, col=colByChr(cs.gr, colors = "brewer.set1"))
kpPlotRainfall(kp, data = cs.gr, col=colByChr(cs.gr, colors = "brewer.set1"),ymax=8)
#kpAxis(kp, ymax = 7, tick.pos = 1:7)
kpAxis(kp, ymax=8, tick.pos = 1:8)
kpAddLabels(kp, labels = c("Distance between markers (log10)"), srt=90, pos=1, label.margin = 0.07)


## Create a combined Rainfall and density plot with custom plotting parameters
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20
pp$leftmargin <- 0.08
kp <- plotKaryotype(plot.type=4, genome=equ_genome, chromosomes=c(1:31,"X","Y"),
                    labels.plotter = NULL, plot.params = pp)
kpAddChromosomeNames(kp, srt=90, cex=0.7)
kpPlotRainfall(kp, data = cs.gr, col=colByChr(cs.gr, colors = "brewer.set1"),ymax=8, r0=0, r1=0.7)
kpAxis(kp, ymax=8, tick.pos = 1:8, r0=0, r1=0.7, cex=0.8)
kpAddLabels(kp, labels = c("Distance between markers (log10)"), srt=90, pos=1, label.margin = 0.07, r0=0, r1=0.7)
kpPlotDensity(kp, data = cs.gr, r0=0.73, r1=1, col=lighter("#889F34"))
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.73, r1=1, cex=0.8)
kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.07, r0=0.71, r1=1)



## Make a Rainfall plot for a zoom-in region
zoom.region <- toGRanges(data.frame("X", 700000, 2071050))
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20
pp$leftmargin <- 0.08
kp <- plotKaryotype(plot.type=4,genome=equ_genome, chromosomes="X", zoom=zoom.region,
                    labels.plotter = NULL, plot.params = pp)
kpAddChromosomeNames(kp, srt=90, cex=0.7)
kpPlotRainfall(kp, data = cs.gr, ymax=8)
kpAxis(kp, ymax=8, tick.pos = 1:8, cex=0.8)
kpAddLabels(kp, labels = c("Distance between markers (log10)"), srt=90, pos=1, label.margin = 0.07)