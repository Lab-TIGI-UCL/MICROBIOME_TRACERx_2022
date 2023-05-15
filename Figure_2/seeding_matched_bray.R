# Microbial landscape of NSCLC 
# Author - Krupa Thakkar 
# Date - 19th October 2022 

# loading packages 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(viridis)
library(vegan)

# Loading data 
absolute_counts <- read.csv("~/Files/normalised_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Tx_data/taxmat.rds") %>% as.data.frame()

region_level <- read.delim("~/seedingRegionInfo.txt", header = TRUE) %>% as.data.frame()

############
# Step 1: Large dataframe based on primary and nonprimary within and between patients  
############
meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(!Biological == "Normal")
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL")
#meta_data <- meta_data %>% filter(!Patient == "U_LTX682")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

plist <- c("Tumour", "Lymph_Node")
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
full_region_samples <- region_level$PrimaryRegion

for (i in 1:nrow(full_region_mat)) {
  dat <- full_region_mat[i, ] %>% as.data.frame()
  
  if (dat$Type_1 == "Primary") {
    
    primary_sample <- dat$Sample_1 %>% as.character()
    nonprimary_sample <- dat$Sample_2 %>% as.character()
    
    if (primary_sample %in% full_region_samples) {
      
      primary_sample <- primary_sample %>% as.character()
      
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
      
    } else {
      full_region_mat[i, "seeding"] <- "Not_Present"
    }}
  
  else {
    primary_sample <- dat$Sample_2 %>% as.character()
    nonprimary_sample <- dat$Sample_1 %>% as.character()
    
    if (primary_sample %in% full_region_samples) {
      
      primary_sample <- primary_sample %>% as.character()
      
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
    } else {
      full_region_mat[i, "seeding"] <- "Not_Present"
    }}
  print(paste0(i, "_complete"))}


full_region_mat <- full_region_mat %>% as.data.frame()

for (i in 1:nrow(full_region_mat)) {
  dat <- full_region_mat[i, ] %>% as.data.frame()
  sample_1 <- dat$Sample_1
  sample_2 <- dat$Sample_2
  
  batch_1 <- meta_data %>% filter(Sample == sample_1)
  batch_1 <- batch_1$Batch 
  batch_2 <- meta_data %>% filter(Sample == sample_2)
  batch_2 <- batch_2$Batch 
  
  full_region_mat[i, "batch_1"] <- batch_1
  full_region_mat[i, "batch_2"] <- batch_2
  
  if (batch_1 == batch_2) {
    full_region_mat[i, "batch_same"] <- TRUE 
  }  else {
    full_region_mat[i, "batch_same"] <- FALSE
  }
}

absolute_counts[absolute_counts < 10] <- 0
absolute_counts <- absolute_counts %>% t() %>% as.data.frame()
braycurtis = vegdist(absolute_counts, "bray")
bray<-as.matrix (braycurtis)
bray <- bray %>% as.data.frame()
bray$samples <- rownames(bray)

bray <- bray %>% melt()
bray <- bray %>% as.data.frame()

test <- full_region_mat %>% filter(!seeding == "Not_Present")

test$bray <- NA

for (i in 1:nrow(test)) {
  dat <- test[i, ] %>% as.data.frame()
  sample_1 <- dat$Sample_1 %>% as.character()
  sample_2 <- dat$Sample_2 %>% as.character()
  
  tmp_bray <- bray %>% filter(samples == sample_1)
  tmp_bray <- tmp_bray %>% filter(variable == sample_2) %>% as.data.frame()
  val <- tmp_bray$value
  
  test[i, "bray"] <- val
  
}

order <- c("TRUE", "FALSE")
test$seeding <- factor(test$seeding, levels = order)

p <- ggplot(test, aes(x = seeding, y = bray, fill = seeding)) + geom_violin(width = 1) + geom_boxplot(width = 0.1, color = "white", alpha = 0.2) + theme_classic() + stat_compare_means() + ylab("Bray Curtis scores") + xlab("Primary seeding metastasis")
p + scale_fill_manual(values = c("deepskyblue4", "darkorange", "deeppink")) + theme(axis.line = element_line(colour = "black"), 
                                                                                    axis.title.x=element_text(size = rel(2)), 
                                                                                    axis.title.y=element_text(size = rel(2)),
                                                                                    legend.text=element_text(size=rel(2)), 
                                                                                    axis.text.y = element_text(face="bold",size=14),
                                                                                    axis.text.x = element_text(face="bold",size=14))

tmp_new <- test %>% select(seeding, bray, batch_same)
order <- c("TRUE", "FALSE")
tmp_new$seeding <- factor(tmp_new$seeding, levels = order)
tmp_new$batch_same <- factor(tmp_new$batch_same, levels = order)

fit = lm(bray ~ batch_same + seeding, data = tmp_new)
summary(fit)

nonprimary_bray <- read.csv("~/Files/nonprimary_bray.csv") %>% as.data.frame()
rownames(nonprimary_bray) <- nonprimary_bray[, 1] 
nonprimary_bray[, 1] <- NULL
nonprimary_bray <- nonprimary_bray %>% as.data.frame()

test <- nonprimary_bray %>% na.omit()

test$bray <- NA

for (i in 1:nrow(test)) {
  dat <- test[i, ] %>% as.data.frame()
  sample_1 <- dat$Sample_1 %>% as.character()
  sample_2 <- dat$Sample_2 %>% as.character()
  
  tmp_bray <- bray %>% filter(samples == sample_1)
  tmp_bray <- tmp_bray %>% filter(variable == sample_2) %>% as.data.frame()
  val <- tmp_bray$value
  
  test[i, "bray"] <- val
  
}

order <- c("INTRA", "EXTRA")
test$location_new <- factor(test$location_new, levels = order)

p <- ggplot(test, aes(x = location_new, y = bray, fill = location_new)) + geom_violin(width = 1) + geom_boxplot(width = 0.1, color = "white", alpha = 0.2) + theme_classic() + stat_compare_means() + ylab("Bray Curtis scores") + xlab("Site of mets")
p + scale_fill_manual(values = c("deepskyblue4", "darkorange", "deeppink")) + theme(axis.line = element_line(colour = "black"), 
                                                                                    axis.title.x=element_text(size = rel(2)), 
                                                                                    axis.title.y=element_text(size = rel(2)),
                                                                                    legend.text=element_text(size=rel(2)), 
                                                                                    axis.text.y = element_text(face="bold",size=14),
                                                                                    axis.text.x = element_text(face="bold",size=14))

tmp_new <- test %>% select(location_new, bray, batch_same)
order <- c("TRUE", "FALSE")
tmp_new$batch_same <- factor(tmp_new$batch_same, levels = order)

fit = lm(bray ~ batch_same + location_new, data = tmp_new)
summary(fit)



# matched vs unmatched primary non primary analysis 
to_keep <- c("UM_Primary_NonPrimary","M_Primary_NonPrimary" )
small_region_mat <- tmp %>% filter(Full_Name %in% to_keep)

absolute_counts[absolute_counts < 10] <- 0
absolute_counts <- absolute_counts %>% t() %>% as.data.frame()
braycurtis = vegdist(absolute_counts, "bray")
bray<-as.matrix (braycurtis)
bray <- bray %>% as.data.frame()
bray$samples <- rownames(bray)

bray <- bray %>% melt()
bray <- bray %>% as.data.frame()

small_region_mat$bray <- NA

for (i in 1:nrow(small_region_mat)) {
  dat <- small_region_mat[i, ] %>% as.data.frame()
  sample_1 <- dat$Sample_1 %>% as.character()
  sample_2 <- dat$Sample_2 %>% as.character()
  
  tmp_bray <- bray %>% filter(samples == sample_1)
  tmp_bray <- tmp_bray %>% filter(variable == sample_2) %>% as.data.frame()
  val <- tmp_bray$value
  
  small_region_mat[i, "bray"] <- val
  
}

for (i in 1:nrow(small_region_mat)) {
  dat <- small_region_mat[i, ] %>% as.data.frame()
  sample_1 <- dat$Sample_1
  sample_2 <- dat$Sample_2
  
  batch_1 <- meta_data %>% filter(Sample == sample_1)
  batch_1 <- batch_1$Batch 
  batch_2 <- meta_data %>% filter(Sample == sample_2)
  batch_2 <- batch_2$Batch 
  
  small_region_mat[i, "batch_1"] <- batch_1
  small_region_mat[i, "batch_2"] <- batch_2
  
  if (batch_1 == batch_2) {
    small_region_mat[i, "batch_same"] <- TRUE 
  }  else {
    small_region_mat[i, "batch_same"] <- FALSE
  }
}

order <- c("TRUE", "FALSE")
test$seeding <- factor(test$seeding, levels = order)

p <- ggplot(small_region_mat, aes(x = Full_Name, y = bray, fill = Full_Name)) + geom_violin(width = 1) + geom_boxplot(width = 0.1, color = "white", alpha = 0.2) + theme_classic() + stat_compare_means() + ylab("Bray Curtis scores") + xlab("Primary and nonprimary")
p + scale_fill_manual(values = c("deepskyblue4", "darkorange", "deeppink")) + theme(axis.line = element_line(colour = "black"), 
                                                                                    axis.title.x=element_text(size = rel(2)), 
                                                                                    axis.title.y=element_text(size = rel(2)),
                                                                                    legend.text=element_text(size=rel(2)), 
                                                                                    axis.text.y = element_text(face="bold",size=14),
                                                                                    axis.text.x = element_text(face="bold",size=14))

