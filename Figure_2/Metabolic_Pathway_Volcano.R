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

##### Step 3 : Assignment of the pathways offline to make the volcano plot #####
pathway_mat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Misc/kegg_pathways_method3.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
pathway_mat$Pathways_2[pathway_mat$Pathways_2==""] <- NA

pathways_hl <- table(pathway_mat$Pathways_2) %>% as.data.frame()
pathways_hl <- pathways_hl %>% filter(Freq > 2)
pathways_hl <- pathways_hl$Var1

pathway_mat$final_pathway <- NA

for (i in 1:nrow(pathway_mat)) {
  
  dat <- pathway_mat[i, ] %>% as.data.frame()
  pathway <- dat$Pathways_2 %>% as.character()
  
  if(pathway %in% pathways_hl) {
    
    pathway_mat[i, "final_pathway"] <- pathway
    
  } else {
    
    pathway_mat[i, "final_pathway"] <- NA
    
  }
}

P_38 <- c("#F7001C", "#2EFB16", "#222AF5", "#F8CACB", "#FF00E2", "#00FAFC", "#F6E60D", "#16660D", "#F58826", "#2A4D7F", "#C41669", "#EB9FFE", 
          "#5DFAAE", "#763B16", "#C400F7", "#0DABFE", "#D1E79F", "#7D2668", "#C1DBFF","#6D22AD","#E4B168","#6D7671","#A3DE16","#FF8FC8", 
          "#D4382E", "#9BDFCA", "#EA817F", "#AA00A6", "#FF0075", "#766E16", "#0094A5", "#A29BFA", "#FB22B5", "#B28BB6", "#22A31C", "#26D4FD", "#0058C4","#9F65FD")

p <- ggplot(data=pathway_mat, aes(x=foldchange, y=-1*log10(pvalue), col=final_pathway)) + geom_point() + theme_minimal()
p + theme_classic() + scale_color_manual(values = P_38) + xlab("Fold Change") + ylab("-log10(p_value)") + theme(
  axis.text.x=element_text(size = rel(2)),
  axis.ticks.x=element_blank(), 
  panel.background = element_blank(), 
  axis.text.y = element_text(size = rel(2)), 
  legend.text=element_text(size=rel(1)), 
  axis.title.y = element_text(size = rel(2)), 
  axis.title.x = element_text(size = rel(2)))
