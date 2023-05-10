# Microbiome of TRACERx NSCLC - Tumour enriched bacteria bar plot and ubiquitous box plot 
# Author - Krupa Thakkar 
# Date - 1st February 2023 

# Loading libraries 
library(tidyverse)
library(ggplot2)
library(reshape2)
library(EnhancedVolcano)
library(DESeq2)
library(nlme)
library(compositions)
library(readr)

# Loading data 
absolute_counts <- read.csv("~/Processed/voom_snm_normalised_new.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Meta_data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

##########
# Step 1 : Tumour enriched bacteria bar plot
##########
meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% select(samples_to_keep) %>% as.data.frame()


clonality_mat <- matrix(nrow = 23, ncol = 4) %>% as.data.frame()
otu_of_interest <- c("Acidovorax", "Aerococcus", "Brevundimonas", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", 
                     "Exiguobacterium", "Flaviflexus", "Fusobacterium", "Hymenobacter", "Klebsiella", "Lachnoanaerobaculum", 
                     "Lactococcus", "Parapusillimonas", "Pseudomonas", "Rothia", "Streptococcus", 
                     "Tepidiphilus", "Thauera", "unidentified_Enterobacteriaceae", "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae")
rownames(clonality_mat) <- otu_of_interest
colnames(clonality_mat) <- c("ubiquitous", "not_ubiquitous", "absent", "single_region")

patient_list <- unique(meta_data$Patient)
absolute_counts <- absolute_counts %>% t() %>% as.data.frame()
absolute_counts$patient <- NA

for (i in rownames(absolute_counts)) {
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Sample == i)
  patient <- tmp_meta$Patient %>% as.character()
  absolute_counts[i, "patient"] <- patient
  
}

for (i in rownames(clonality_mat)) {
  i <- i %>% as.character()
  tmp_mat <- absolute_counts %>% select(i, "patient")
  
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
clonality_mat <- clonality_mat[order(clonality_mat$ubiquitous), ]

order_of_genus <- rownames(clonality_mat)

clonality_mat <- clonality_mat %>% melt()
clonality_mat <- clonality_mat %>% as.data.frame()
clonality_mat$genus <- factor(clonality_mat$genus,
                              levels = order_of_genus)
values <- c("absent", "single_region","not_ubiquitous", "ubiquitous")
clonality_mat$variable <- factor(clonality_mat$variable,
                                 levels = values)


t<- ggplot(clonality_mat, aes(fill=variable, y=value, x=genus)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("darkorange","darkviolet", "deeppink", "deepskyblue4")) + coord_flip() + ylab("Proportion of Patients") + xlab("Tumour Associated Genus") + theme_classic()
t + theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          panel.background = element_blank(), 
          axis.title.x = element_text(size = rel(2)), 
          axis.title.y = element_text(size = rel(2)), 
          legend.text=element_text(size=rel(1)), 
          axis.text.y = element_text(face="bold",size=14))
t  

##### Are bacteria that are present ubiquitously more likely to be tumour enriched ##### 
microbes <- deseq2_results$name
clonality_mat <- matrix(nrow = 109, ncol = 4) %>% as.data.frame()
rownames(clonality_mat) <- microbes
colnames(clonality_mat) <- c("ubiquitous", "not_ubiquitous", "absent", "single_region")

for (i in rownames(clonality_mat)) {
  i <- i %>% as.character()
  tmp_mat <- absolute_counts %>% select(i, "patient")
  
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
clonality_mat$name <- rownames(clonality_mat)

full_dat <- merge(clonality_mat, deseq2_results, by = "name")
full_dat <- full_dat %>% select(-"tumour_enriched")
vals <- c("UP", "DOWN")
full_dat <- full_dat %>% filter(diffexpressed %in% vals)

clonality_dat$diffexpressed <- factor(clonality_dat$diffexpressed,
                                      levels = vals)

p <- ggplot(full_dat, aes(x = diffexpressed, y = ubiquitous, fill = diffexpressed)) + geom_boxplot()
p + stat_compare_means() 
p + theme_classic() + theme(axis.title.x = element_text(size = rel(1.5)), 
                            axis.title.y = element_text(size = rel(1.5)), 
                            legend.text=element_text(size=rel(1.5)), 
                            axis.text.y = element_text(face="bold",size=14),
                            axis.text.x = element_text(face="bold",size=14)) + xlab("Tumour enriched") + ylab("Proportion of Patients") + scale_fill_manual(values = c("lightpink","orchid"))
