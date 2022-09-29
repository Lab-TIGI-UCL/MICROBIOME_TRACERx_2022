# Microbial landscape of NSCLC - Smoking associated microbes and pathways 
# Author - Krupa Thakkar 
# Date - 22nd September 2022 

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(DESeq2)
library(nlme)
library(compositions)
library(readr)

# Loading data
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/new_samples_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/tax_mat.csv") %>% as.data.frame()

##### Step 1: Data formatting 
meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

tmp <- absolute_counts %>% as.data.frame()
tmp[tmp < 10] <- 0
tmp[tmp > 9] <- 1
tmp$total <- rowSums(tmp)
#tmp <- tmp %>% filter(total > 20)
genus_to_keep <- rownames(tmp)

absolute_counts <- subset(absolute_counts, rownames(absolute_counts) %in% genus_to_keep)

##### Step 2 : Adjusting the meta data for smoking vs never smokers (similar group sizes) #####
meta_data[meta_data == "Never Smoked"] <- "Never_Smoker"
meta_data[meta_data == "Recent Ex-Smoker"] <- "Ever_Smoker"
meta_data[meta_data == "Current Smoker"] <- "Ever_Smoker"

ns <- meta_data %>% filter(Smoking_Status == "Never_Smoker")
s <- meta_data %>% filter(Smoking_Status == "Ever_Smoker")
new_test <- rbind(ns, s)
new_test <- new_test %>% as.data.frame()

absolute_ns <- absolute_counts %>% select(ns$Sample)
tmp_ns <- absolute_ns
tmp_ns[tmp_ns < 10] <- 0
tmp_ns[tmp_ns > 9] <- 1
tmp_ns <- tmp_ns %>% as.data.frame()
tmp_ns$total <- rowSums(tmp_ns)
tmp_ns <- tmp_ns %>% filter(total > 1)
ns_otu <- rownames(tmp_ns)
absolute_ns <- subset(absolute_ns, rownames(absolute_ns) %in% ns_otu)
absolute_ns$otu <- rownames(absolute_ns)

absolute_s <- absolute_counts %>% select(s$Sample)
tmp_s <- absolute_s
tmp_s[tmp_s < 10] <- 0
tmp_s[tmp_s > 9] <- 1
tmp_s <- tmp_s %>% as.data.frame()
tmp_s$total <- rowSums(tmp_s)
tmp_s <- tmp_s %>% filter(total > 8)
s_otu <- rownames(tmp_s)
absolute_s <- subset(absolute_s, rownames(absolute_s) %in% s_otu)
absolute_s$otu <- rownames(absolute_s)

absolute_dat <- merge(absolute_ns, absolute_s, by = "otu", all = TRUE) %>% as.data.frame()
rownames(absolute_dat) <- absolute_dat[, 1]
absolute_dat[, 1] <- NULL
absolute_dat <- absolute_dat %>% as.data.frame()
absolute_dat[is.na(absolute_dat)] <- 0
absolute_dat <- absolute_dat %>% as.data.frame()

smoking_counts <- absolute_dat %>% as.data.frame()
smoking_counts <- smoking_counts + 1

##### Step 2 : Differential abundance  for smokers and never smokers #####
new_test$Smoking_Status <- as.factor(new_test$Smoking_Status)

dds <- DESeqDataSetFromMatrix(countData = smoking_counts, 
                              colData = new_test, 
                              design = ~ Smoking_Status)
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)

head(results(dds, tidy=TRUE))
deseq2_pvals <- res@listData[["padj"]] %>% as.data.frame()
deseq2_foldchange <- res@listData[["log2FoldChange"]] %>% as.data.frame()
deseq2_results <- cbind(deseq2_pvals, deseq2_foldchange)
colnames(deseq2_results) <- c("pvalue", "foldchange")
rownames(deseq2_results) <- rownames(results)
deseq2_results$probename <- rownames(deseq2_results)
deseq2_results$diffexpressed <- "NO"
deseq2_results$diffexpressed[deseq2_results$foldchange > 0.5 & deseq2_results$pvalue < 0.05] <- "UP"
deseq2_results$diffexpressed[deseq2_results$foldchange < -0.5 & deseq2_results$pvalue < 0.05] <- "DOWN"
deseq2_results$name <- rownames(smoking_counts)

p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p + theme_classic() + scale_color_manual(values = c("darkorange", "deeppink", "deepskyblue4")) + xlab("Fold Change") 

