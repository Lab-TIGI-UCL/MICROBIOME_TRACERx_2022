# Microbial landscape of NSCLC - Microbial enrichment in hot and cold immune regions 
# Author - Krupa Thakkar 
# Date - 15th September 2022 

# Loading libraries 
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(EnhancedVolcano)
library(DESeq2)
library(nlme)
library(compositions)
library(readr)

# Loading data 
absolute_counts <- read.csv("~/Files/decontaminated_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

load("/Data/Additional/tx100_immune_KL.RData")

regional <- read.csv("~/DAdditional/regional_immmune.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(regional) <- regional[, 1] 
regional[, 1] <- NULL
regional <- regional %>% as.data.frame()

##### Step 1: Formatting the data to keep only the regions with immune and micorbial data and then subsetting the meta_data to inlcude this (Also in this analysis the intermediate will be moved to the cold) ####
meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample

absolute_counts[absolute_counts < 10] <- 0
absolute_counts <- absolute_counts%>% select(samples_to_keep)
regional <- regional %>% filter(Sample %in% samples_to_keep)
regional <- regional %>% select(-name)
#regional[regional == "intermediate"] <- "cold"
regional <- regional %>% filter(!immuneClass_2 == "intermediate")


samples_to_keep <- regional$Sample
absolute_counts <- absolute_counts%>% select(samples_to_keep)
absolute_counts <- absolute_counts %>% as.data.frame()

tmp_normal <- absolute_counts
tmp_normal[tmp_normal < 10] <- 0
tmp_normal[tmp_normal > 9] <- 1
tmp_normal <- tmp_normal %>% as.data.frame()
tmp_normal$total <- rowSums(tmp_normal)
tmp_normal <- tmp_normal %>% filter(total > 4)
normal_otu <- rownames(tmp_normal)
absolute_counts <- subset(absolute_counts, rownames(absolute_counts) %in% normal_otu)
absolute_counts <- absolute_counts %>% as.data.frame()

rownames(regional) <- regional$Sample
regional$immuneClass_2 <- as.factor(regional$immuneClass_2)
absolute_counts <- round(absolute_counts)
absolute_counts <- absolute_counts + 1

regional$batches <- NA

for (i in rownames(regional)) {
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Sample == i)
  batch <- tmp_meta$Batch
  
  regional[i, "batches"] <- batch
}

##### Step 2 : Differential abundance analysis #####
dds <- DESeqDataSetFromMatrix(countData = absolute_counts, 
                              colData = regional, 
                              design = ~ batches + immuneClass_2)
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
deseq2_results$name <- rownames(absolute_counts)


##### Step 3 : Plotting data #####
p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p + theme_classic() + scale_color_manual(values = c("darkturquoise","green", "darkviolet"))

to_keep <- c("UP", "DOWN")
deseq2_results <- deseq2_results %>% filter(diffexpressed %in% to_keep)
to_keep <- c("Escherichia-Shigella", "Atopostipes", 'Pantoea', "Streptomyces", "Campylobacter", "Clostridia_UCG-014", "unidentified_Micrococcaceae", 
             "Tissierella", "Peptoniphilus", "Alistipes"
             )
deseq2_results <- deseq2_results %>% filter(name %in% to_keep)
deseq2_results <- deseq2_results %>% as.data.frame()
deseq2_results$name <- factor(deseq2_results$name,
                              levels = to_keep)

color <- ifelse(deseq2_results$foldchange < 0, "darkorange", "darkviolet")
p <- ggplot(deseq2_results, aes(x = reorder(name, foldchange), y = foldchange)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,      # Background color
           color = "white") + # Remove the legend
  xlab("Group") +
  ylab("Fold Change ") + coord_flip() + theme_classic()
p + theme(axis.title.x=element_text(size = rel(2)),
          axis.text.x=element_text(size = rel(2)),
          axis.ticks.x=element_blank(), 
          axis.title.y = element_blank(),
          panel.background = element_blank(), 
          axis.text.y = element_text(size = rel(2)), 
          legend.text=element_text(size=rel(2)))




