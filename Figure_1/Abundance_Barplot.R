# Microbial landscape of NSCLC - Abundance and prevalence plots 
# Author : Krupa Thakkar 
# Date - 15th September 2022 

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Loading data
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()
absolute_counts[absolute_counts < 10] <- 0

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

##### Step 1 : Abundance plot #####
meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% select(samples_to_keep) %>% as.data.frame()
phylum_counts <- absolute_counts %>% as.data.frame()
phylum_counts <- phylum_counts %>% na.omit()
phylum_counts$phylum <- NA

taxmat_genus <- taxmat$genus

for (i in rownames(phylum_counts)) {
  
  i <- i %>% as.character()
  
  if (i %in% taxmat_genus) {
    taxmat_tmp <- taxmat %>% filter(genus == i)
    taxmat_tmp <- taxmat_tmp[1, ] %>% as.data.frame()
    phylum <- taxmat_tmp$phylum
    phylum_counts[i, "phylum"] <- phylum
  } else {
    
    family <- strsplit(i, split = "unidentified_")
    family <- family[[1]][2] %>% as.character()
    
    taxmat_tmp <- taxmat %>% filter(family == family)
    taxmat_tmp <- taxmat_tmp[1, ] %>% as.data.frame()
    phylum <- taxmat_tmp$phylum
    phylum_counts[i, "phylum"] <- phylum
    
  }
}

phylum_counts <- phylum_counts %>% as.data.frame()

to_keep <- c("Proteobacteria", "Firmicutes", "Actinobacteriota", "Campilobacterota", "Bacteroidota", "Spirochaetota", "Methylomirabilota")

for (i in rownames(phylum_counts)) {
  
  phylum <- phylum_counts[i, ] %>% as.data.frame()
  phylum <- phylum$phylum %>% as.character()
  
  
  if (phylum %in% to_keep) {
    
    phylum_counts[i, "phylum"] <- phylum
  } else {
    
    phylum_counts[i, "phylum"] <- "Other"
    
  }
}

phylum_counts <- phylum_counts %>% group_by(phylum) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
rownames(phylum_counts) <- phylum_counts$phylum
phylum_counts <- phylum_counts %>% select(-phylum)
phylum_counts <- phylum_counts %>% as.data.frame()
phylum_counts <- phylum_counts %>% t() %>% as.data.frame()
phylum_counts$patient <- NA

for (i in rownames(phylum_counts)) {
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Sample == i) %>% as.data.frame()
  patient <- tmp_meta$Patient %>% as.character()
  
  phylum_counts[i, "patient"] <- patient
}

phylum_counts <- phylum_counts %>% as.data.frame()
phylum_counts <- phylum_counts %>% group_by(patient) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
phylum_counts <- phylum_counts[order(-phylum_counts$Proteobacteria), ]

patient_order <- phylum_counts$patient

phylum_melt <- phylum_counts %>% melt() %>% as.data.frame()

phylum_melt_tmp <- phylum_counts
rownames(phylum_melt_tmp) <- phylum_melt_tmp$patient
phylum_melt_tmp <- phylum_melt_tmp %>% select(-patient) %>% as.data.frame()
phylum_melt_tmp$total <- rowSums(phylum_melt_tmp)
phylum_melt_tmp <- (phylum_melt_tmp/phylum_melt_tmp[,9])*100
phylum_melt_tmp <- phylum_melt_tmp %>% select(-total) %>% as.data.frame()
phylum_melt_tmp <- phylum_melt_tmp[order(-phylum_melt_tmp$Proteobacteria), ]
phylum_melt_tmp$patient <- rownames(phylum_melt_tmp)
patient_order <- phylum_melt_tmp$patient

phylum_melt <- phylum_melt_tmp %>% melt() %>% as.data.frame()


order_of_phylum <- c("Other", "Methylomirabilota", "Spirochaetota", "Bacteroidota", "Campilobacterota", "Actinobacteriota", "Firmicutes", "Proteobacteria")
phylum_melt$patient <- factor(phylum_melt$patient,
                              levels = patient_order)
phylum_melt$variable <- factor(phylum_melt$variable,
                               levels = order_of_phylum)

colours <- c( "yellow", "darkorange", "darkolivegreen", "darkseagreen3", "darkslategray2", "darkviolet",   "deeppink4","deepskyblue4"  )



p <- ggplot(phylum_melt, aes(fill = variable, y = value, x = patient)) + geom_bar(position="fill", stat="identity")
p + scale_fill_manual(values = colours) + theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(), 
                                                panel.background = element_blank())
