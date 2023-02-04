library(tidyverse)
library(ggpubr)
setwd("/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-snpdensity/pruned_dataset/femaleRef/snpden")

# MALE
## load in files
myFiles <- list.files(pattern="Male*")
## build backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
## load other files
for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}
SNPdensity.males <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

# FEMALE
## load in files
myFiles <- list.files(pattern="Female*")
## build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
## loop over other files
for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}
SNPdensity.females <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")


# Combine male and female
SNPdensity <- full_join(SNPdensity.males, SNPdensity.females)

# make new tibble, change type for base, and replace NAs
SNPdensity.rows <- SNPdensity
SNPdensity.rows$base <- as.numeric(as.character(SNPdensity.rows$base))
SNPdensity.rows <- SNPdensity.rows %>% replace_na(list(Males = 0, Females = 0)) %>% replace(is.na(.), 0)

# use a subsetter to get the male and female average densities 
male_n <- ncol(SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ])
males_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

female_n <- ncol(SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ])
females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity <- full_join(males_true,females_true) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity$mean_MvF_dif > 0) #15831 pruned=18594
table(true_SNPdensity$mean_MvF_dif < 0) #22225 pruned=22001
table(true_SNPdensity$mean_MvF_dif == 0) #4388 pruned=1849


# initial run of the permutation:
perm <- sample(3:length(SNPdensity.rows))
SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)
males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
true_SNPdensity_position_perm <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
perm_backbone <- true_SNPdensity_position_perm %>% select(scaf,base,mean_MvF_dif) %>% rename(p1=mean_MvF_dif)

## run loops of perm test
for(i in 2:1000){
  perm <- sample(3:length(SNPdensity.rows))
  SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)
  males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
  females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
  column_ID <- paste0("p",i)
  perm_upgrade <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females) %>% select(scaf,base,mean_MvF_dif) %>% dplyr::rename(!!column_ID := "mean_MvF_dif")
  perm_backbone <- full_join(perm_backbone,perm_upgrade)
}
perm_with_true <- full_join(true_SNPdensity %>% select(scaf,base,mean_MvF_dif),perm_backbone)
perm_with_true_p <- perm_with_true %>% mutate(Pvalue = case_when(mean_MvF_dif > 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1), mean_MvF_dif < 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1), mean_MvF_dif == 0 ~ 1))

# chage working directory to save things and access index reference genome
setwd("/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-snpdensity/pruned_dataset/femaleRef")

# write to files
write_tsv(perm_with_true_p, "SNPdensity_perm_with_true_p.femaleRef.sex.txt")
write_tsv(perm_with_true_p %>% select(scaf,base,mean_MvF_dif,Pvalue), "SNPdensity_perm_with_true_p_select.femaleRef.sex.txt")

# want a look at proportion of scaffold - find the outlier
genome_index <- read_tsv("CHRR_integrated_F.fa.fai",col_names = F) %>% rename(scaf=X1,length=X2)
perm_with_true_p001_snp <- perm_with_true_p %>% filter(Pvalue <= 0.001) 
count_snp_p001_scaffolds <- perm_with_true_p001_snp %>% select(scaf) %>% count(scaf)
proportion_count_snp_p001_scaffolds <- left_join(count_snp_p001_scaffolds,genome_index) %>% mutate(proportion = (n*10000)/length)
proportion_count_snp_p001_scaffolds %>% filter(length > 10000000) %>% ggplot(aes(x=scaf,y=proportion)) + geom_col() + theme(axis.text.x = element_text(angle = 90))
ggsave("sex_scaf_proportion_p001.pdf")

proportion_count_snp_p001_scaffolds %>% filter(length > 1000000) %>% ggplot(aes(x=scaf,y=proportion)) + geom_col() + theme(axis.text.x = element_text(angle = 90))
ggsave("sex_scaf_proportion_p001_small_1000000.pdf")
