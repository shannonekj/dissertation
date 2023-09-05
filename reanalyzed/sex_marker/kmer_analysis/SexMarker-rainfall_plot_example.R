#!/bin/R
# Taken from https://github.com/drtamermansour/Equine_parentage_testing/blob/main/scripts/rainfall.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
#BiocManager::install("Hmisc")


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