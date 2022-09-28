# Microbial landscape of NSCLC - HE, Spatial analysis for immune cell interactions 
# Author - Krupa Thakkar 
# Date - 28th September 2022

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(vegan)

# Loading data 
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

load("/Users/krupathakkar/Documents/Projects/Microbial_TRACERx_P1/Data/Additional/tx100_immune_KL.RData")

regional <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/Additional/regional_immmune.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(regional) <- regional[, 1] 
regional[, 1] <- NULL
regional <- regional %>% as.data.frame()

tx <- tx %>% as.data.frame()

#### Step 1 :  Diversity, unique gneus and microbial load metrics ##### 
absolute_counts[absolute_counts < 10] <- 0

meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

absolute_diversity <- absolute_counts %>% t() %>% as.data.frame()
absolute_diversity <- diversity(absolute_diversity, index = "shannon") %>% as.data.frame()
colnames(absolute_diversity) <- "diversity"
absolute_diversity$names <- rownames(absolute_diversity)

meta_data$diversity <- NA
meta_data$unique_genus <- NA
meta_data$microbial_load <- NA
  
for (i in samples_to_keep) {
    
    i <- i %>% as.character()
    
    tmp_dat <- absolute_counts %>% select(i)
    
    tmp_dat[tmp_dat == 0] <- NA
    tmp_dat <- tmp_dat %>% na.omit()
    
    unique_genus <- nrow(tmp_dat)
    meta_data[i, "unique_genus"] <- unique_genus
    
    tmp_dat <- tmp_dat %>% t() %>% as.data.frame()
    tmp_dat$total <- rowSums(tmp_dat)
    microbial_load <- tmp_dat$total
    meta_data[i, "microbial_load"] <- microbial_load
    
    diversity <- absolute_diversity %>% filter(names == i)
    diversity <- diversity$diversity
    meta_data[i, "diversity"] <- diversity
    
  }

regional <- regional %>% filter(Sample %in% samples_to_keep)

meta_data <- merge(meta_data, regional, by = "Sample", all = TRUE)
meta_data <- meta_data[!is.na(meta_data$file_name), ]
meta_data <- meta_data %>% as.data.frame()

meta_data$load_class <- NA
for (i in 1:nrow(meta_data)) {
  tmp_dat <- meta_data[i, ] %>% as.data.frame()
  load <- tmp_dat$microbial_load 
  
  if (load < 291) {
    meta_data[i, "load_class"] <- "low"
  } else {
    meta_data[i, "load_class"] <- "high"
  }
}

comparisons <- c("hot", "cold")
cbPalette2 <- c("deepskyblue", "deepskyblue4")
p <- ggplot(meta_data, aes(x = load_class, y = fi)) + geom_boxplot(aes(fill = load_class), outlier.shape = NA) + geom_jitter(color="black", size=1, alpha=0.5) + theme_classic()
p + stat_compare_means() + scale_fill_manual(values = c("deeppink", "deepskyblue4")) + xlab("Microbial Load Class")

