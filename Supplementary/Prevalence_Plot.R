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

##### Step 2 : Prevalence plots + ubiquitous plots####
genus_prevalence <- absolute_counts
binary_prevalence <- genus_prevalence
binary_prevalence[binary_prevalence < 10] <- 0
binary_prevalence[binary_prevalence > 9] <- 1
binary_prevalence$total <- rowSums(binary_prevalence)
binary_prevalence <- binary_prevalence[order(-binary_prevalence$total), ]
binary_prevalence <- binary_prevalence[1:30, ]
binary_prevalence <- binary_prevalence %>% as.data.frame()
binary_prevalence <- binary_prevalence %>% select(-total)
genus_to_keep <- rownames(binary_prevalence)

genus_prevalence$genus <- rownames(genus_prevalence)
genus_prevalence <- genus_prevalence %>% filter(genus %in% genus_to_keep )
genus_prevalence <- genus_prevalence %>% select(-genus) %>% as.data.frame()
genus_prevalence <- genus_prevalence %>% t() %>% as.data.frame()
genus_prevalence$patient <- NA

for (i in rownames(genus_prevalence)) {
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Sample == i)
  patient <- tmp_meta$Patient %>% as.character()
  genus_prevalence[i, "patient"] <- patient
  
}

patient_list <- genus_prevalence$patient %>% unique()

clonality_mat <- matrix(nrow = nrow(binary_prevalence), ncol = 4) %>% as.data.frame()
rownames(clonality_mat) <- rownames(binary_prevalence)
colnames(clonality_mat) <- c("ubiquitous", "not_ubiquitous", "single_region", "absent")

for (i in rownames(clonality_mat)) {
  i <- i %>% as.character()
  tmp_mat <- genus_prevalence %>% select(i, "patient")
  
  patient_counts <- matrix(nrow = 220, ncol = 1) %>% as.data.frame()
  rownames(patient_counts) <- patient_list
  colnames(patient_counts) <- c("value")
  
  for (j in patient_list) {
    tmp_2 <- tmp_mat %>% filter(patient == j)
    tmp_3 <- tmp_2 %>% select(-patient)
    tmp_3[tmp_3 < 10] <- 0
    tmp_3[tmp_3 > 9] <- 1
    tmp_3 <- tmp_3 %>% t() %>% as.data.frame()
    total <- rowSums(tmp_3)
    
    if (total > 0) {
      total_samples <- nrow(tmp_2)
      
      if (total_samples < 2) {
        patient_counts[j, "value"] <- "single_region"
      } else {
        new_value <- total/total_samples
        
        if (new_value == 1) {
          patient_counts[j, "value"] <- "ubiquitous"
        } else {
          patient_counts[j, "value"] <- "not_ubiquitous"
        }
      }} else {
        patient_counts[j, "value"] <- "absent"
      }
    
  }
  patient_counts <- patient_counts %>% as.data.frame()
  shared <- patient_counts %>% filter(value == "ubiquitous")
  shared <- (nrow(shared)/220)*100
  clonality_mat[i, "ubiquitous"] <- shared 
  
  patient_counts <- patient_counts %>% as.data.frame()
  single <- patient_counts %>% filter(value == "single_region")
  single <- (nrow(single)/220)*100
  clonality_mat[i, "single_region"] <- single
  
  not_shared <- patient_counts %>% filter(value == "not_ubiquitous")
  not_shared <- (nrow(not_shared)/220)*100
  clonality_mat[i, "not_ubiquitous"] <- not_shared 
  
  absent <- patient_counts %>% filter(value == "absent")
  absent <- (nrow(absent)/220)*100
  clonality_mat[i, "absent"] <- absent 
  
}

clonality_mat <- clonality_mat %>% as.data.frame()
clonality_mat$genus <- rownames(clonality_mat)
# to melt the clonality mat 

binary_prevalence <- binary_prevalence %>% as.data.frame()
binary_prevalence <- binary_prevalence %>% t() %>% as.data.frame()
binary_prevalence$patient <- NA
for (i in rownames(binary_prevalence)) {
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Sample == i)
  patient <- tmp_meta$Patient %>% as.character()
  binary_prevalence[i, "patient"] <- patient
  
}

binary_prevalence <- binary_prevalence %>% group_by(patient) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
rownames(binary_prevalence) <- binary_prevalence[, 1] 
binary_prevalence[, 1] <- NULL
binary_prevalence <- binary_prevalence %>% as.data.frame()
binary_prevalence[binary_prevalence > 0] <- 1
binary_prevalence <- binary_prevalence %>% t() %>% as.data.frame()
binary_prevalence$total <- rowSums(binary_prevalence)
binary_prevalence$genus <- rownames(binary_prevalence)
binary_prevalence <- binary_prevalence %>% select(total, genus)
binary_prevalence$prop <- "prop"
binary_prevalence <- binary_prevalence %>% as.data.frame()
binary_prevalence <- binary_prevalence[order(-binary_prevalence$total), ]
order_of_genus <- binary_prevalence$genus

binary_prevalence$genus <- factor(binary_prevalence$genus,
                                  levels = order_of_genus)


q <- ggplot(binary_prevalence, aes(x=genus, y=prop, size = total, color = prop)) +
  geom_point(alpha=0.7) + theme_classic() + 
  theme(legend.position="bottom") +
  ylab("Prevalence") +
  xlab("Genus") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="bold",size=14), 
        axis.text.y = element_text(face="bold",size=14))  + scale_color_manual(values = "darkgreen")
q

clonality_mat <- clonality_mat %>% melt()
clonality_mat <- clonality_mat %>% as.data.frame()
clonality_mat$genus <- factor(clonality_mat$genus,
                                  levels = order_of_genus)
values <- c("absent", "single_region", "not_ubiquitous", "ubiquitous")
clonality_mat$variable <- factor(clonality_mat$variable,
                              levels = values)


t<- ggplot(clonality_mat, aes(fill=variable, y=value, x=genus)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("darkorange","darkviolet", "deeppink", "deepskyblue4")) 
t + theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          panel.background = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(2)), 
          legend.position = "none", 
          axis.text.y = element_text(face="bold",size=14))
t 
