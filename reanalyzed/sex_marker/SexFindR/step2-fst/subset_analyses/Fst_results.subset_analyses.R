# Fst Analysis - Both Reference Genome
## with subset of individuals that aren't outliers

### DIRECTIONS:
#       - open in the directory containing *maleRef directories corresponding to Fst files from subset of individuals


##################
###   R WORK   ###
##################
library(tidyverse)
setwd("/group/millermrgrp2/shannon/projects/dissertation/reanalyzed/sex_marker/SexFindR/step2-fst/subset_analyses")

##### --- Female Reference Genome --- #####
Fst_F <- read_tsv("femaleRef/femaleRef_no_outliers_pruned_fst.weir.fst") %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0)) %>% rename(scaf = CHROM, base=POS) %>% mutate(base=as.numeric(base)) %>% filter(WEIR_AND_COCKERHAM_FST >=0) %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0))
scaffold_lengths_F <- read_tsv("../../step0-variantcalling/femaleRef/CHRR_integrated.fa.fai", col_names = c("scaf","length")) %>% filter(grepl("lg",scaf)) 
cutoff1_F <- round(nrow(Fst_F)*0.05)
Fst_F_sort_cut1 <- head(Fst_F %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff1_F) 
min(Fst_F_sort_cut1$`WEIR_AND_COCKERHAM_FST`) # lowest value there is 0.118254  (5% cutoff)
min_five_cutoff_F <- min(Fst_F_sort_cut1$`WEIR_AND_COCKERHAM_FST`)
cutoff2_F <- round(nrow(Fst_F)*0.01)
Fst_F_sort_cut2 <- head(Fst_F %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff2_F) 
min(Fst_F_sort_cut2$`WEIR_AND_COCKERHAM_FST`) # lowest value there is 0.193863 (1% cutoff)
min_one_cutoff_F <- min(Fst_F_sort_cut2$`WEIR_AND_COCKERHAM_FST`)
table(Fst_F_sort_cut2$scaf)
quantile(Fst_F$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)
mean(Fst_F$WEIR_AND_COCKERHAM_FST) # 0.03776015
Fst_female <- left_join(scaffold_lengths_F,Fst_F) %>% filter(WEIR_AND_COCKERHAM_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR_F <- length(unique(Fst_female$scaf))
Fst_female$BPcum <- 0
sf <- 0
nbp_F <- c()
for (i in unique(Fst_female$scaf)){
  nbp_F[i] <- max(Fst_female[Fst_female$scaf == i,]$base)
  Fst_female[Fst_female$scaf == i,"BPcum"] <- Fst_female[Fst_female$scaf == i,"base"] + sf
  sf <- sf + nbp_F[i]
}
axis.set_F <- Fst_female %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
Fst_female %>% ggplot(aes(x=BPcum, y=WEIR_AND_COCKERHAM_FST, color=as.factor(scaf))) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis.set_F$scaf, breaks = axis.set_F$center) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="CHROMOSOME", color="") +
  guides(color = guide_legend(ncol = 1, byrow = F)) +
  geom_hline(yintercept = min_five_cutoff_F, linetype="dotted", color = "black", size=0.75) +
  geom_hline(yintercept = min_one_cutoff_F, linetype="dotted", color = "black", size=1.5)  +
  geom_vline(xintercept = as.integer(sf), linetype="dotted", color="red", size=1) +
  labs(title="Fst for delta smelt (5% and 1% cutoff & FemaleSexSDR as dotted lines)")
ggsave("femaleRef_no_outliers_Fst_manhattan_with95_99_noMissing.pdf")




##### --- Male Reference Genome --- #####
Fst_M <- read_tsv("maleRef/maleRef_no_outliers_pruned.weir.fst") %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0)) %>% rename(scaf = CHROM, base=POS) %>% mutate(base=as.numeric(base)) %>% filter(WEIR_AND_COCKERHAM_FST >=0) %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0))
scaffold_lengths_M <- read_tsv("../../step0-variantcalling/maleRef/CHRR_integrated.fa.fai", col_names = c("scaf","length")) %>% filter(grepl("lg",scaf)) 
cutoff1_M <- round(nrow(Fst_M)*0.05)
Fst_M_sort_cut1 <- head(Fst_M %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff1_M) 
min(Fst_M_sort_cut1$`WEIR_AND_COCKERHAM_FST`) # lowest value there is 0.120286 (5% cutoff)
min_five_cutoff_M <- min(Fst_M_sort_cut1$`WEIR_AND_COCKERHAM_FST`)
cutoff2_M <- round(nrow(Fst_M)*0.01)
Fst_M_sort_cut2 <- head(Fst_M %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff2_M) 
min(Fst_M_sort_cut2$`WEIR_AND_COCKERHAM_FST`) # lowest value there is 0.186137 (1% cutoff)
min_one_cutoff_M <- min(Fst_M_sort_cut2$`WEIR_AND_COCKERHAM_FST`)
table(Fst_M_sort_cut2$scaf)
quantile(Fst_M$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)
mean(Fst_M$WEIR_AND_COCKERHAM_FST) # 0.03747152
Fst_male <- left_join(scaffold_lengths_M, Fst_M) %>% filter(WEIR_AND_COCKERHAM_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR_M <- length(unique(Fst_male$scaf))
Fst_male$BPcum <- 0
sm <- 0
nbp_M <- c()
for (i in unique(Fst_male$scaf)){
  nbp_M[i] <- max(Fst_male[Fst_male$scaf == i,]$base)
  Fst_male[Fst_male$scaf == i,"BPcum"] <- Fst_male[Fst_male$scaf == i,"base"] + sm
  sm <- sm + nbp_M[i]
}
axis.set_M <- Fst_male %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
Fst_male %>% ggplot(aes(x=BPcum, y=WEIR_AND_COCKERHAM_FST, color=as.factor(scaf))) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis.set_M$scaf, breaks = axis.set_M$center) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="CHROMOSOME", color="") +
  guides(color = guide_legend(ncol = 1, byrow = F)) +
  geom_hline(yintercept = min_five_cutoff_M, linetype="dotted", color = "black", size=0.75) +
  geom_hline(yintercept = min_one_cutoff_M, linetype="dotted", color = "black", size=1.5)  +
  geom_vline(xintercept = as.integer(sm), linetype="dotted", color="red",size=1) +
  labs(title="Fst for delta smelt (5% and 1% cutoff & MaleSexSDR as dotted lines)")
ggsave("maleRef_no_outliers_Fst_manhattan_with95_99_noMissing.pdf")


