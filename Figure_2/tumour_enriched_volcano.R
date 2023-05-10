# Microbiome of TRACERx NSCLC - Tumour enriched bacteria volcano plot 
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
absolute_counts <- read.csv("~/Files/16S_data/Processed/decontaminated_final.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Files/Meta_data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

##########
# Step 1 - Creating two dataframes for tumour vs normal
########## 
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL")
meta_data <- meta_data%>% filter(!Batch == "2.2")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

otu_tables <- rownames(absolute_counts) %>% as.data.frame()
colnames(otu_tables) <- "full_name"
otu_tables$value <- rownames(otu_tables)
otu_tables$value<- sub("^", "otu_", otu_tables$value)

normal_meta <- meta_data %>% filter(Biological == "Normal")
normal_samples <- normal_meta$Sample
normal_patients <- normal_meta$Patient
absolute_normal <- absolute_counts %>% select(normal_samples)

tumour_meta  <- meta_data %>% filter(Biological == "Tumour")
tumour_meta <- tumour_meta %>% filter(Patient %in% normal_patients)
tumour_samples <- tumour_meta$Sample
tumour_patients <- tumour_meta$Patient
absolute_tumour <- absolute_counts %>% select(tumour_samples)

tmp_normal <- absolute_normal
tmp_normal[tmp_normal < 10] <- 0
tmp_normal[tmp_normal > 9] <- 1
tmp_normal <- tmp_normal %>% as.data.frame()
tmp_normal$total <- rowSums(tmp_normal)
tmp_normal <- tmp_normal %>% filter(total > 3)
normal_otu <- rownames(tmp_normal)
absolute_normal <- subset(absolute_normal, rownames(absolute_normal) %in% normal_otu)
absolute_normal$otu <- rownames(absolute_normal)

tmp_tumour_paired <- absolute_tumour
tmp_tumour_paired[tmp_tumour_paired < 10] <- 0
tmp_tumour_paired[tmp_tumour_paired > 9] <- 1
tmp_tumour_paired <- tmp_tumour_paired %>% as.data.frame()
tmp_tumour_paired$total <- rowSums(tmp_tumour_paired)
tmp_tumour_paired <- tmp_tumour_paired %>% filter(total > 7)
tumour_otu_paired <- rownames(tmp_tumour_paired)
absolute_tumour <- subset(absolute_tumour, rownames(absolute_tumour) %in% tumour_otu_paired)
absolute_tumour <- absolute_tumour %>% t() %>% as.data.frame()
absolute_tumour$patient <- NA

for (i in rownames(absolute_tumour)) {
  i <- i %>% as.character()
  tmp <- meta_data %>% filter(Sample == i) %>% as.data.frame()
  patient <- tmp$Patient %>% as.character()
  absolute_tumour[i, "patient"] <- patient
}

absolute_tumour <- absolute_tumour %>% group_by(patient) %>%
  dplyr::summarise(across(everything(), mean))
absolute_tumour <- absolute_tumour%>% as.data.frame()
rownames(absolute_tumour) <- absolute_tumour$patient
absolute_tumour <- absolute_tumour %>% dplyr::select(-patient)
absolute_tumour <- absolute_tumour %>% t() %>% as.data.frame()
colnames(absolute_tumour)<-paste(colnames(absolute_tumour),"T",sep="-")
absolute_tumour$otu <- rownames(absolute_tumour)

absolute_normal <- absolute_normal %>% dplyr::select(-otu)
absolute_normal <- absolute_normal%>% t() %>% as.data.frame()
absolute_normal$patient <- NA

for (i in rownames(absolute_normal)) {
  i <- i %>% as.character()
  tmp <- meta_data %>% filter(Sample == i) %>% as.data.frame()
  patient <- tmp$Patient %>% as.character()
  absolute_normal[i, "patient"] <- patient
}

absolute_normal <- absolute_normal %>% group_by(patient) %>%
  dplyr::summarise(across(everything(), mean))
absolute_normal <- absolute_normal %>% as.data.frame()
rownames(absolute_normal) <- absolute_normal$patient
absolute_normal <- absolute_normal %>% dplyr::select(-patient)
absolute_normal <- absolute_normal %>% t() %>% as.data.frame()
colnames(absolute_normal)<-paste(colnames(absolute_normal),"N",sep="-")
absolute_normal$otu <- rownames(absolute_normal)

######## 
# Step 2 - Merging both the tumour and normal dataframes 
#######
absolute_dat_paired <- merge(absolute_normal, absolute_tumour, by = "otu", all = TRUE) %>% as.data.frame()
rownames(absolute_dat_paired) <- absolute_dat_paired[, 1]
absolute_dat_paired[, 1] <- NULL
absolute_dat_paired <- absolute_dat_paired %>% as.data.frame()
absolute_dat_paired[is.na(absolute_dat_paired)] <- 0
absolute_dat_paired <- absolute_dat_paired %>% as.data.frame()

absolute_paired_meta <- colnames(absolute_dat_paired) %>% as.data.frame()
colnames(absolute_paired_meta) <- "sample"
rownames(absolute_paired_meta) <- absolute_paired_meta$sample
absolute_paired_meta$patient <- NA
for (i in rownames(absolute_paired_meta)) {
  tmp <- i %>% as.character()
  tmp <- strsplit(tmp, split = "-")[[1]][1]
  absolute_paired_meta[i, "patient"] <- tmp
}
absolute_paired_meta <- absolute_paired_meta %>% as.data.frame()
absolute_paired_meta$biological <- NA
for (i in rownames(absolute_paired_meta)) {
  tmp <- i %>% as.character()
  tmp <- strsplit(tmp, split = "-")[[1]][2]
  absolute_paired_meta[i, "biological"] <- tmp
}
absolute_paired_meta <- absolute_paired_meta %>% as.data.frame()
absolute_paired_meta$biological <- as.factor(absolute_paired_meta$biological)

absolute_dat_paired <- round(absolute_dat_paired)
absolute_dat_paired <- absolute_dat_paired + 1 # Adding a psuedocount of 1 

##########
# Step 3 : DESEQ2 for differential abudnance analysis 
##########
dds <- DESeqDataSetFromMatrix(countData = absolute_dat_paired, 
                              colData = absolute_paired_meta, 
                              design = ~ patient + biological)
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
deseq2_results$name <- rownames(absolute_dat_paired)
deseq2_results$tumour_enriched <- NA

p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p + theme_classic()

genus_of_interest <- c("Acidovorax", "Aerococcus", "Brevundimonas", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", 
                       "Exiguobacterium", "Flaviflexus", "Fusobacterium", "Hymenobacter", "Klebsiella", "Lachnoanaerobaculum", 
                       "Lactococcus", "Parapusillimonas", "Providencia", "Pseudomonas", "Rothia", "Streptococcus", 
                       "Tepidiphilus", "Thauera", "unidentified_Enterobacteriaceae", "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae")
for (i in 1:nrow(deseq2_results)) {
  
  tmp <- deseq2_results[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  
  if (name %in% genus_of_interest) {
    
    deseq2_results[i, "tumour_enriched"] <- name
    
  } else {
    
    deseq2_results[i, "tumour_enriched"] <- NA
  }
}
deseq2_results <- deseq2_results %>% as.data.frame()

P_38 <- c("#F7001C", "#2EFB16", "#222AF5", "#F8CACB", "#FF00E2", "#00FAFC", "#F6E60D", "#16660D", "#F58826", "#2A4D7F", "#C41669", "#EB9FFE", 
          "#5DFAAE", "#763B16", "#C400F7", "#0DABFE", "#D1E79F", "#7D2668", "#C1DBFF","#6D22AD","#E4B168","#6D7671","#A3DE16","#FF8FC8", 
          "#D4382E", "#9BDFCA", "#EA817F")


p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed, label = tumour_enriched)) + geom_point(size = 3, alpha = 0.6) + theme_classic() + theme(axis.title.x = element_text(size = rel(2)), 
                                                                                                                                                                                  axis.title.y = element_text(size = rel(2)), 
                                                                                                                                                                                  legend.position="none", 
                                                                                                                                                                                  axis.text.y = element_text(face="bold",size=14),
                                                                                                                                                                                  axis.text.x = element_text(face="bold",size=14)) + scale_color_manual(values = c("deeppink", "darkorange", "deepskyblue3"))
p + geom_label_repel(aes(label=ifelse(diffexpressed == "UP", as.character(name), '')), box.padding= 0.35, point.padding = 0.9,segment.color = 'grey50', max.overlaps = 20, color = "black") + xlab("Fold Change (Tumour vs Normal)") 




