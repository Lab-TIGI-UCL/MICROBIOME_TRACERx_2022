# Microbial landscape of NSCLC - Seeding vs non-seeding bacteria 
# Author - Krupa Thakkar 
# Date - 2nd March 2023 

# Loading libraries 
library(tidyverse)
library(dplyr)
library(reshape2)
library(EnhancedVolcano)
library(DESeq2)
library(nlme)
library(compositions)
library(readr)

# Loading data 
matrix <- read.csv("~/Files/decontaminated_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(matrix) <- matrix[, 1] 
matrix[, 1] <- NULL
matrix <- matrix %>% as.data.frame()
counts <- matrix

seeding <- read.csv("~/Files/seeding_bray.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(seeding) <- seeding[, 1] 
seeding[, 1] <- NULL
seeding <- seeding %>% as.data.frame()

# Step 1 : Reformatting data 
# Samples to test are LTX469 R2/R3 for seeding and R4 for non-seeding 
primary_samples_1 <- seeding %>% filter(Type_1 == "Primary")
primary_samples_1 <- primary_samples_1$Sample_1
primary_samples_2 <- seeding %>% filter(Type_2 == "Primary")
primary_samples_2 <- primary_samples_2$Sample_2

full_list <- c(primary_samples_1, primary_samples_2) %>% unique()

new_meta <- matrix(nrow = 42, ncol = 1) %>% as.data.frame()
rownames(new_meta) <- full_list
colnames(new_meta) <- "seeding"
new_meta$seeding <- c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "TRUE", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "FALSE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "FALSE", "TRUE", "TRUE", "TRUE")
new_meta <- new_meta %>% as.data.frame()
new_meta$seeding <- as.factor(as.character(new_meta$seeding))
new_meta$sample <- rownames(new_meta)
samples_to_keep <- rownames(new_meta)

counts <- counts %>% select(samples_to_keep)

test <- counts %>% as.data.frame()
test[test < 10] <- 0
test[test > 9] <- 1
test$total <- rowSums(test)
test <- test %>% filter(total > 2)
genus_to_keep <- rownames(test)

counts <- subset(counts, rownames(counts) %in% genus_to_keep)


counts <- counts + 1 

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = new_meta, 
                              design = ~ seeding)
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
deseq2_results$name <- rownames(counts)
deseq2_results$tumour_enriched <- NA

p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed)) + geom_point() + geom_point(size = 3, alpha = 0.6) + theme_classic() + theme(axis.title.x = element_text(size = rel(2)), 
                                                                                                                                                                        axis.title.y = element_text(size = rel(2)), 
                                                                                                                                                                        legend.position="none", 
                                                                                                                                                                        axis.text.y = element_text(face="bold",size=14),
                                                                                                                                                                        axis.text.x = element_text(face="bold",size=14)) + scale_color_manual(values = c("deeppink", "darkorange", "deepskyblue3"))
p + theme_classic()
















