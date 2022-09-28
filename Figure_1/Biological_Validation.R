# Microbial landscape of NSCLC - Pneumonia biological validation plots 
# Author - Krupa Thakkar 
# Date - 26th September 2022 

# Loading libraries and packages 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Loading data 
matrix <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(matrix) <- matrix[, 1] 
matrix[, 1] <- NULL
matrix <- matrix %>% as.data.frame()
counts <- matrix

meta_data <- read.csv("~/Documents/QIIME2_pipeline/Tx_data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

#full_dat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/association_meta_new.csv", header = TRUE) %>% as.data.frame()

old_meta <- read.csv("~/Documents/Project 2 -  Microbial landscape of NSCLC/Data/sample_meta_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame() 

##### Step 1 : Reformatting data and creating a large meta data to include all of the metrics of interest #####
batches <- c("2.2")
meta_data <- meta_data%>% filter(!Batch %in% batches)
samples_to_keep <- c("Tumour")
meta_data <- meta_data %>% filter(Biological %in% samples_to_keep)
samples_to_keep <- meta_data$Sample

counts <- counts %>% dplyr::select(samples_to_keep)
counts[counts < 10] <- 0

full_dat <- colnames(counts) %>% as.data.frame()
colnames(full_dat) <- "samples"
rownames(full_dat) <- full_dat$samples
full_dat$unique_genus <- NA
full_dat$microbial_load <- NA
full_dat$staphylococcus <- NA
full_dat$pseudomonas <- NA
full_dat$streptococcus <- NA
full_dat$klebsiella <- NA
full_dat$total_score <- NA
full_dat$emphysema <- NA
full_dat$copd <- NA

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp <- counts %>% select(i) %>% as.data.frame()
  
  unique_genus <- tmp %>% as.data.frame()
  unique_genus[unique_genus == 0] <- NA
  unique_genus <- unique_genus %>% na.omit()
  unique_genus <- nrow(unique_genus) 
  full_dat[i, "unique_genus"] <- unique_genus
  
  microbial_load <- tmp %>% as.data.frame()
  microbial_load[microbial_load == 0] <- NA
  microbial_load <- microbial_load %>% na.omit()
  microbial_load <- colSums(microbial_load) 
  full_dat[i, "microbial_load"] <- microbial_load
  
  bacteria_tmp <- tmp %>% t() %>% as.data.frame()
  bacteria_tmp <- bacteria_tmp %>% select("Staphylococcus", "Streptococcus", "Pseudomonas", "Klebsiella")
  staphy <- bacteria_tmp$Staphylococcus
  full_dat[i, "staphylococcus"] <- staphy
  strepto <- bacteria_tmp$Streptococcus
  full_dat[i, "streptococcus"] <- strepto
  kleb <- bacteria_tmp$Klebsiella
  full_dat[i, "klebsiella"] <- kleb
  
  bacteria_tmp[bacteria_tmp < 10] <- NA
  bacteria_tmp <- bacteria_tmp %>% t() %>% as.data.frame()
  bacteria_tmp <- bacteria_tmp %>% na.omit() %>% as.data.frame()
  bacteria_tmp <- nrow(bacteria_tmp) 
  full_dat[i, "total_score"] <- bacteria_tmp
  
}
full_dat <- full_dat %>% as.data.frame()
full_dat$pathology <- NA

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp_meta <- meta_data %>% filter(Sample == i)
  pathology <- tmp_meta$PathTissueFind 
  full_dat[i, "pathology"] <- pathology
}

full_dat$FEV1 <- NA

for (i in 1:nrow(full_dat)) {
  
  tmp <- full_dat[i, ] %>% as.data.frame()
  sample <- tmp$samples %>% as.character()
  
  fev1 <- old_meta[old_meta$Sample %in% sample, ] 
  fev1 <- fev1[1, ] %>% as.data.frame()
  fev1 <- fev1$FEV1_group
  
  full_dat[i, "FEV1"] <- fev1

}

full_dat <- full_dat %>% as.data.frame()

##### Step 2 : Plotting genus level associations ##### 
full_dat[is.na(full_dat)] <- FALSE
order <- c("TRUE", "FALSE")
full_dat$pneumonia <- factor(full_dat$pneumonia, 
                             level = order)
full_dat$total_score <- as.numeric(as.character(full_dat$total_score))
full_dat$emphysema <- factor(full_dat$emphysema, 
                             level = order)

tmp <- table(full_dat$emphysema, full_dat$total_score) %>% as.data.frame()
tmp$Var1 <- factor(tmp$Var1, level = order)

p <- ggplot(full_dat, aes(x = FEV1, y = streptococcus, fill = FEV1)) + geom_boxplot()
p +  theme_classic() + stat_compare_means() + scale_fill_brewer(palette = "BuPu") + xlab("FEV1 group") + ylab("Streptococcus Counts (Normalised)")

p <- ggplot(tmp, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity", position = "fill")
p +  theme_classic()+ scale_fill_manual(values = c("lightpink","orchid","deeppink","deepskyblue4")) + xlab("Emphysema Diagnosis") + ylab("Total pathogenic genus")
