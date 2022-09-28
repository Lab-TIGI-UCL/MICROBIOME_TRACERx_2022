# Microbial Landscape of NSCLC - Seeding level analysis (primary /non-primary)
# Author - Krupa Thakkar 
# Date - 28th September 2022 

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(viridis)

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

region_level <- read.delim("~/Documents/QIIME2_pipeline/seedingRegionInfo.txt", header = TRUE) %>% as.data.frame()

###### Step 1: Large dataframe based on primary and nonprimary within and between patients  #####
meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(!Biological == "Normal")
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

plist <- c("Tumour", "Lymph_Node")
#plist <- c("Tumour") # To test if the interaction is through tumours or lymph nodes 
primary_list <- meta_data %>% filter(Biological %in% plist)
primary_list <- primary_list$Sample
nonprimary_list <- meta_data %>% filter(!Biological %in% plist)
nonprimary_list <- nonprimary_list$Sample

# Full data frame 
total_samples <- colnames(absolute_counts)
tmp <- t(combn(total_samples, 2))
tmp <- tmp %>% as.data.frame()
colnames(tmp) <- c("Sample_1", "Sample_2")
tmp$Patient_1 <- NA
tmp$Type_1 <- NA
tmp$Patient_2 <- NA
tmp$Type_2 <- NA

primary <- primary_list %>% as.data.frame()
colnames(primary) <- "Sample"
primary$Type <- "Primary"

nonprimary <- nonprimary_list %>% as.data.frame()
colnames(nonprimary) <- "Sample"
nonprimary$Type <- "NonPrimary"

sample_info <- rbind(primary, nonprimary)
rownames(sample_info) <- sample_info$Sample
sample_info$Patient <- NA
for (i in rownames(sample_info)) {
  i <- i %>% as.character()
  dat <- meta_data %>% filter(Sample == i)
  p <- dat$Patient %>% as.character()
  sample_info[i, "Patient"] <- p
}


for (i in 1:nrow(tmp)) {
  dat <- tmp[i, ] %>% as.data.frame()
  s1 <- dat$Sample_1 %>% as.character()
  test_1 <- sample_info %>% filter(Sample == s1)
  p_1 <- test_1$Patient %>% as.character()
  t_1 <- test_1$Type %>% as.character()
  tmp[i, "Patient_1"] <- p_1
  tmp[i, "Type_1"] <- t_1
  
  s2 <- dat$Sample_2 %>% as.character()
  test_2 <- sample_info %>% filter(Sample == s2)
  p_2 <- test_2$Patient %>% as.character()
  t_2 <- test_2$Type %>% as.character()
  tmp[i, "Patient_2"] <- p_2
  tmp[i, "Type_2"] <- t_2
}

tmp$same <- NA
for (i in 1:nrow(tmp)) {
  r1 <- dat$Sample_1 %>% as.character()
  r2 <- dat$Sample_2 %>% as.character()
  if (r1 == r2) {
    tmp[i, "same"] <- "YES"
  }
  else {
    tmp[i, "same"] <- "NO"
  }
}

tmp$total <- NA
for (i in 1:nrow(tmp)) {
  dat <- tmp[i, ] %>% as.data.frame()
  p_1 <- dat$Patient_1 %>% as.character()
  p_2 <- dat$Patient_2 %>% as.character()
  
  if (p_1 == p_2) {
    tmp[i, "total"] <- "Matched"
  }
  else {
    tmp[i, "total"] <- "Unmatched"
  } 
}

full_region_mat <- tmp %>% as.data.frame()
full_region_mat <- full_region_mat %>% filter(total == "Matched")

my_cols <- c("total", "Type_1", "Type_2")
full_region_mat$Full_Name <- do.call(paste, c(full_region_mat[my_cols], sep = "_"))
full_region_mat[full_region_mat == "Matched_NonPrimary_Primary"] <- "M_Primary_NonPrimary"
full_region_mat[full_region_mat == "Matched_NonPrimary_NonPrimary"] <- "M_NonPrimary_NonPrimary"
full_region_mat[full_region_mat == "Matched_Primary_NonPrimary"] <- "M_Primary_NonPrimary"
full_region_mat[full_region_mat == "Matched_Primary_Primary"] <- "M_Primary_Primary"

full_region_mat <- full_region_mat %>% filter(Full_Name == "M_Primary_NonPrimary")

full_region_mat$seeding <- NA

seeding_samples <- region_level %>% filter(Metastasizing == TRUE)
seeding_samples <- seeding_samples$PrimaryRegion

for (i in 1:nrow(full_region_mat)) {
  dat <- full_region_mat[i, ] %>% as.data.frame()
  
  if (dat$Type_1 == "Primary") {
    
    primary_sample <- dat$Sample_1 %>% as.character()
    nonprimary_sample <- dat$Sample_2 %>% as.character()
    
    if (primary_sample %in% seeding_samples) {
      
      mets <- region_level %>% filter(PrimaryRegion == primary_sample)
      mets <- mets$Metastasis %>% strsplit(mets, split = ";")
      mets <- mets[[1]]
      
      if (nonprimary_sample %in% mets) {
        full_region_mat[i, "seeding"] <- TRUE
      } else {
        full_region_mat[i, "seeding"] <- FALSE
      }} else {
        full_region_mat[i, "seeding"] <- FALSE
      }
    
  }
  
  else {
    primary_sample <- dat$Sample_2 %>% as.character()
    nonprimary_sample <- dat$Sample_1 %>% as.character()
    
    if (primary_sample %in% seeding_samples) {
      mets <- region_level %>% filter(PrimaryRegion == primary_sample)
      mets <- mets$Metastasis %>% strsplit(mets, split = ";")
      mets <- mets[[1]]
      
      if (nonprimary_sample %in% mets) {
        full_region_mat[i, "seeding"] <- TRUE
      } else {
        full_region_mat[i, "seeding"] <- FALSE
      }} else {
        full_region_mat[i, "seeding"] <- FALSE
      }
  }
}

for (i in 1:nrow(full_region_mat)) {
  
  dat <- full_region_mat[i, ] %>% as.data.frame()
  
  if (dat$Type_1 == "Primary") {
    
    primary_sample <- dat$Sample_1 %>% as.character()
    nonprimary_sample <- dat$Sample_2 %>% as.character()
    
    if (primary_sample %in% seeding_samples) {
      
      mets <- region_level %>% filter(PrimaryRegion == primary_sample)
      mets <- mets$Metastasis %>% strsplit(mets, split = ";")
      mets <- mets[[1]] 
      
      if (nonprimary_sample %in% mets) {
        full_region_mat[i, "seeding"] <- TRUE
      } else {
        full_region_mat[i, "seeding"] <- FALSE
      }
    } else {
      full_region_mat[i, "seeding"] <- FALSE
    }
  } else {
    
    primary_sample <- dat$Sample_2 %>% as.character()
    nonprimary_sample <- dat$Sample_1 %>% as.character()
    
    if (primary_sample %in% seeding_samples) {
      
      mets <- region_level %>% filter(PrimaryRegion == primary_sample)
      mets <- mets$Metastasis %>% strsplit(mets, split = ";")
      mets <- mets[[1]] 
      
      if (nonprimary_sample %in% mets) {
        full_region_mat[i, "seeding"] <- TRUE
      } else {
        full_region_mat[i, "seeding"] <- FALSE
      }
    } else {
      full_region_mat[i, "seeding"] <- FALSE
    }
  }
  
  print(paste0(i, "_complete")) }


full_region_mat <- full_region_mat %>% as.data.frame()

for (i in 1:nrow(full_region_mat)) {
  dat <- full_region_mat[i, ]
  
  s_1 <- dat$Sample_1 %>% as.character()
  s_2 <- dat$Sample_2 %>% as.character()
  
  s1_matrix <- absolute_counts %>% select(s_1)
  s1_matrix <- s1_matrix %>% as.data.frame()
  s1_matrix[s1_matrix < 10] <- NA
  s1_matrix <- na.omit(s1_matrix)
  s1_matrix$genus <- rownames(s1_matrix)
  s1_tmp <- rownames(s1_matrix)
  
  s2_matrix <- absolute_counts %>% select(s_2)
  s2_matrix <- s2_matrix %>% as.data.frame()
  s2_matrix[s2_matrix < 10] <- NA
  s2_matrix <- na.omit(s2_matrix)
  s2_matrix$genus <- rownames(s2_matrix) 
  s2_tmp <- rownames(s2_matrix)
  
  if (dat$Type_1 == "Primary") {
    s2_matrix <- s2_matrix %>% filter(genus %in% s1_tmp)
    s2_tmp <- rownames(s2_matrix)
    j_score <- length(intersect(s1_tmp,s2_tmp))/length(union(s1_tmp,s2_tmp))
    full_region_mat[i, "jaccard_score"] <- j_score
  } else {
    s1_matrix <- s1_matrix %>% filter(genus %in% s2_tmp)
    s1_tmp <- rownames(s1_matrix)
    j_score <- length(intersect(s1_tmp,s2_tmp))/length(union(s1_tmp,s2_tmp))
    full_region_mat[i, "jaccard_score"] <- j_score
  }
}

##### Step 2 : Plotting the data #####
full_region_mat <- full_region_mat %>% as.data.frame()
full_region_mat <- full_region_mat %>% na.omit()
order <- c("TRUE","FALSE")
full_region_mat$seeding <- factor(full_region_mat$seeding, 
                                  levels = order)

p <- ggplot(full_region_mat, aes(x = seeding, y = jaccard_score, fill = seeding))  +geom_violin(width=1) + geom_boxplot(width=0.1, color="white", alpha=0.2) + stat_compare_means()
p + scale_fill_brewer(palette = "Paired") + xlab("Region seeding metastasis") + ylab("Microbiome Genus Similarity")
p + theme(panel.background = element_blank())

p <- ggplot(full_region_mat, aes(x = seeding, y = jaccard_score, fill = seeding)) + geom_violin(width = 1) + geom_boxplot(width = 0.1, color = "white", alpha = 0.2) + theme_classic()
p + scale_fill_manual(values = c("deepskyblue4", "darkorange")) + theme( axis.text.x = element_text(size = 10), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.title.y=element_blank())
