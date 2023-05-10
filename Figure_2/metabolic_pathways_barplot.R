# Microbial landscape of NSCLC - Metabolic pathways 16S 
# Author - Krupa Thakkar 
# Date - 5th April 2023 

# Loading libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)

# Loading data 
metabolic_pathway_full <- read.csv("~/MISC/metabolic_pathways_da.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
metabolic_pathway_full[,8] <- NULL
metabolic_pathway_full[,1] <- NULL
up_pathways <- metabolic_pathway_full %>% filter(diffexpressed == "UP")
up_pathways_all <- up_pathways$Pathway %>% unique()
down_pathways <- metabolic_pathway_full %>% filter(diffexpressed == "DOWN")
down_pathways_all <- down_pathways$Pathway %>% unique()

up_specific <- setdiff(up_pathways_all, down_pathways_all)
up_specific <- up_specific %>% as.data.frame()
colnames(up_specific) <- "pathway"
up_specific$direction <- "UP"
down_specific <- setdiff(down_pathways_all, up_pathways_all)
down_specific <- down_specific %>% as.data.frame()
colnames(down_specific) <- "pathway"
down_specific$direction <- "DOWN"

full_dat <- rbind(up_specific, down_specific)
full_dat <- full_dat %>% as.data.frame()
full_dat$fc <- NA

for (i in 1:nrow(full_dat)) {
  tmp <- full_dat[i, ] %>% as.data.frame()
  pathway <- tmp$pathway
  direction <- tmp$direction 
  
  if (direction == "UP") {
    
    tmp_2 <- metabolic_pathway_full[order(-metabolic_pathway_full$foldchange),]
    tmp_2 <- tmp_2 %>% as.data.frame()
    tmp_2 <- tmp_2 %>% filter(Pathway == pathway)
    tmp_2 <- tmp_2[1, ] %>% as.data.frame()
    fc <- tmp_2$foldchange 
    full_dat[i, "fc"] <- fc
    
  } else {
    
    tmp_2 <- metabolic_pathway_full[order(metabolic_pathway_full$foldchange),]
    tmp_2 <- tmp_2 %>% as.data.frame()
    tmp_2 <- tmp_2 %>% filter(Pathway == pathway)
    tmp_2 <- tmp_2[1, ] %>% as.data.frame()
    fc <- tmp_2$foldchange 
    full_dat[i, "fc"] <- fc
    
  }
}

full_dat <- full_dat %>% as.data.frame()
full_dat <- full_dat[order(full_dat$fc),]
rownames(full_dat) <- NULL

p <- ggplot(full_dat, aes( x = reorder(pathway, fc), y = fc, fill = direction)) + geom_bar(stat = "identity")
p + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_fill_manual(values=c("darkorange", "deepskyblue4"))




