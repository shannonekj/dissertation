# Script: plot_PCA-no_outliers.R
# Author: Shannon E.K. Joslin
# Date: 19 January 2023
library(tidyverse)
library(ggpubr)
setwd("/Users/miocene/Desktop/git_repos/dissertation/reanalyzed/sex_marker/SexFindR/step2-gwas/output_no_outliers")

### FEMALE ###
# read in data
pca <- read_table2("./femaleRef_pca.eigenvec", col_names = FALSE)
eigenval <- scan("./femaleRef_pca.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# sort sex
sex <- rep(NA, length(pca$ind))
sex[grep("_F_", pca$ind)] <- "F"
sex[grep("_M_", pca$ind)] <- "M"
# remake data.frame
pca <- as.tibble(data.frame(pca, sex))
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = sex)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("darkred", "darkblue"))
b <- b + coord_equal() 
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
# zoom plot
pca_subset <- pca %>% filter(PC1 < 0.25, PC2 > -0.1, PC2 < 0.25)
c <- ggplot(pca_subset, aes(PC1, PC2, col = sex)) + geom_point(size = 3)
c <- c + scale_colour_manual(values = c("darkred", "darkblue"))
c <- c + coord_equal() 
b+c
# grab outliers
pca_outliers <- pca %>% filter(PC1 > 0.25 | PC2 < -0.1 | PC2 > 0.25)
head(pca_outliers$ind)



### MALE ###
# read in data
pca_m <- read_table2("./maleRef_pca.eigenvec", col_names = FALSE)
eigenval_m <- scan("./maleRef_pca.eigenval")
# sort out the pca data
# remove nuisance column
pca_m <- pca_m[,-1]
# set names
names(pca_m)[1] <- "ind"
names(pca_m)[2:ncol(pca_m)] <- paste0("PC", 1:(ncol(pca_m)-1))
# sort sex
sex <- rep(NA, length(pca_m$ind))
sex[grep("_F_", pca_m$ind)] <- "F"
sex[grep("_M_", pca_m$ind)] <- "M"
# remake data.frame
pca_m <- as.tibble(data.frame(pca_m, sex))
# first convert to percentage variance explained
pve_m <- data.frame(PC = 1:20, pve = eigenval_m/sum(eigenval_m)*100)
# make plot
a_m <- ggplot(pve_m, aes(PC, pve)) + geom_bar(stat = "identity")
a_m + ylab("Percentage variance explained")
# calculate the cumulative sum of the percentage variance explained
cumsum(pve_m$pve)
# plot pca
b_m <- ggplot(pca_m, aes(PC1, PC2, col = sex)) + geom_point(size = 3)
b_m <- b_m + scale_colour_manual(values = c("darkred", "darkblue"))
b_m <- b_m + coord_equal() 
b_m <- b_m + xlab(paste0("PC1 (", signif(pve_m$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_m$pve[2], 3), "%)"))
b_m
# zoom plot
pca_subset_m <- pca_m %>% filter(PC1 > -0.25, PC2 > -0.25)
c_m <- ggplot(pca_subset_m, aes(PC1, PC2, col = sex)) + geom_point(size = 3)
c_m <- c_m + scale_colour_manual(values = c("darkred", "darkblue"))
c_m <- c_m + coord_equal() 
b_m+c_m
# grab outliers
pca_outliers_m <- pca_m %>% filter(PC1 < -0.25 | PC2 < -0.25)
head(pca_outliers_m$ind)



