# Microbial Landscape of NSCLC - Pathway level analysis 
# Author - Krupa Thakkar
# Date - 20th October 2022 

# Loading packages 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(DESeq2)
library(nlme)
library(compositions)
library(readr)
library(gridExtra)

# Loading data 
full_kegg_dat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Misc/metabolic_pathways/full_kegg_dat.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(full_kegg_dat) <- full_kegg_dat[, 1] 
full_kegg_dat[, 1] <- NULL
full_kegg_dat <- full_kegg_dat %>% as.data.frame()

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

##### Step 1 : Building the dataframe to include 5% of tumour and 5% of normal #####
meta_data <- meta_data%>% filter(!Batch == "2.2") 
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL") 
types_to_keep <- c("Tumour", "Normal")
meta_data <- meta_data %>% filter(Biological %in% types_to_keep)
samples_to_keep <- meta_data$Sample 
full_kegg_dat <- full_kegg_dat %>% select(samples_to_keep) 
full_path_dat <- full_path_dat %>% select(samples_to_keep)

normal_data <- meta_data %>% filter(Biological == "Normal") 
normal_samples <- normal_data$Sample 
normal_patients <- normal_data$Patient 
tumour_data <- meta_data %>% filter(Biological == "Tumour") 
tumour_data <- tumour_data %>% filter(Patient %in% normal_patients) 
tumour_samples <- tumour_meta$Sample

to_keep <- c("Tumour", "Normal") 
comb_meta <- rbind(normal_data, tumour_data)
comb_meta <- comb_meta %>% as.data.frame()

full_kegg_dat_normal <- full_kegg_dat %>% select(normal_samples)
tmp_kegg_normal <- full_kegg_dat_normal %>% as.data.frame()
tmp_kegg_normal[tmp_kegg_normal < 100] <- 0 
tmp_kegg_normal[tmp_kegg_normal > 101] <- 1 
tmp_kegg_normal$total <- rowSums(tmp_kegg_normal) 
tmp_kegg_normal <-tmp_kegg_normal %>% filter(total > 3) 
to_keep <- rownames(tmp_kegg_normal)
full_kegg_dat_normal <- subset(full_kegg_dat_normal, rownames(full_kegg_dat_normal) %in% to_keep) 
full_kegg_dat_normal$kegg <- rownames(full_kegg_dat_normal)

full_kegg_dat_tumour <- full_kegg_dat %>% select(tumour_samples)
tmp_kegg_tumour <- full_kegg_dat_tumour %>% as.data.frame()
tmp_kegg_tumour[tmp_kegg_tumour < 100] <- 0 
tmp_kegg_tumour[tmp_kegg_tumour > 101] <- 1 
tmp_kegg_tumour$total <- rowSums(tmp_kegg_tumour) 
tmp_kegg_tumour <- tmp_kegg_tumour %>% filter(total > 7) 
to_keep <- rownames(tmp_kegg_tumour)
full_kegg_dat_tumour <- subset(full_kegg_dat_tumour,rownames(full_kegg_dat_tumour) %in% to_keep) 
full_kegg_dat_tumour$kegg <- rownames(full_kegg_dat_tumour)

kegg_paired <- merge(full_kegg_dat_normal, full_kegg_dat_tumour, by = "kegg",all = TRUE) %>% as.data.frame() 
rownames(kegg_paired) <- kegg_paired[, 1]
kegg_paired[, 1] <- NULL 
kegg_paired <- kegg_paired %>% as.data.frame()
kegg_paired[is.na(kegg_paired)] <- 0 
kegg_paired <- kegg_paired %>% as.data.frame() 
kegg_paired[kegg_paired < 100] <- 0 
kegg_paired <- round(kegg_paired)

tmp_kegg <- kegg_paired + 1 

##### Step 2 : DESEQ2 differential abundance  #####
dds <- DESeqDataSetFromMatrix(countData = tmp_kegg, 
                              colData = comb_meta,
                              design = ~ Batch + Biological) 
dds <- DESeq(dds) 
res <- results(dds,alpha=0.05)

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
deseq2_results$name <- rownames(kegg_paired)

p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p + theme_classic()

# Loading data <- Tumour vs normal pathways assessed using pathways curated from the KEGG API 
tvn_pathways <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Misc/metabolic_pathways/tvn_pathways.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()

# Tumour vs normal specific pathways 
t_pathways <- tvn_pathways$T_Pathways %>% unique()
n_pathways <- tvn_pathways$N_Pathways %>% unique()

t_specific <- setdiff(t_pathways, n_pathways) %>% as.data.frame()
t_specific <- t_specific[1:30, ] %>% as.data.frame()
t_specific$type <- "Tumour"
n_specific <- setdiff(n_pathways, t_pathways) %>% as.data.frame()
n_specific <- n_specific[1:30, ] %>% as.data.frame()
n_specific$type <- "Normal"

full <- rbind(t_specific, n_specific) %>% as.data.frame()

# save to add in the foldchange values 
write.csv(full, "~/Documents/Projects/Microbial_TRACERx_P1/Misc/metabolic_pathways/upregulated_pathways.csv")

full <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Misc/metabolic_pathways/upregulated_pathways.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
color <- ifelse(full$foldchange < 0, "darkblue", "deeppink")
p <- ggplot(full, aes(x = reorder(name, foldchange), y = foldchange)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,      # Background color
           color = "white") + # Remove the legend
  xlab("Group") +
  ylab("Fold Change ")  + theme_classic()
p + theme(axis.title.x=element_text(size = rel(2)),
          axis.text.x=element_text(size = rel(1.5), angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x=element_blank(), 
          axis.title.y = element_blank(),
          panel.background = element_blank(), 
          axis.text.y = element_text(size = rel(2)), 
          legend.text=element_text(size=rel(2)))

