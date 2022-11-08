# Microbial landscape of NSCLC - Nejman et al, Science, 2020 paired NSCLC analysis 
# Author - Krupa Thakkar 
# Date 21st June 2022 

# Loading libraries 
library(tidyverse)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(EnhancedVolcano)
library(readxl)

# Loading data 
taxmat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/tax_mat.csv") %>% as.data.frame()
tx_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/tx_tumour_enriched.csv") %>% as.data.frame()
tx_data <- tx_data$name 

nejman_meta <- read_excel("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Nejman_Science_2020/NEJMAN_META.xlsx") %>% as.data.frame()
rownames(nejman_meta) <- nejman_meta[, 1]
nejman_meta[, 1] <- NULL
nejman_meta <- nejman_meta %>% as.data.frame()

nejman_data <- read_excel("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Nejman_Science_2020/NEJMAN_OTU.xlsx") %>% as.data.frame()
nsclc_meta <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Nejman_Science_2020/Lung_Cancer_META.csv", header = TRUE, check.names = FALSE)
rownames(nsclc_meta) <- nsclc_meta[, 1]
nsclc_meta[, 1] <- NULL
nsclc__meta <- nsclc_meta %>% as.data.frame()

normal_otu <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/normal_otu.csv") %>% as.data.frame()
normal_otu <- normal_otu$x
tumour_otu <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/tumour_otu.csv") %>% as.data.frame()
tumour_otu <- tumour_otu$x

##### Step 1 :  Step 1 : Splitting the meta_data and two dataframes #####
nsclc_meta <- nsclc_meta %>% t() %>% as.data.frame()
samples_to_keep <- rownames(nsclc_meta)

normal_meta <- nsclc_meta %>% filter(Tissue_Type == "Lung (NAT)")
normal_samples <- rownames(normal_meta)
normal_patients <- normal_meta$Pat_ID

tumour_meta <- nsclc_meta %>% filter(Tissue_Type == "Lung (T)")
tumour_meta <- tumour_meta %>% filter(Pat_ID %in% normal_patients)
tumour_patients <- tumour_meta$Pat_ID 

normal_meta <- normal_meta %>% filter(Pat_ID %in% tumour_patients)
tumour_samples <- rownames(tumour_meta)
normal_samples <- rownames(normal_meta)

tumour_meta$Biological <- "Tumour"
tumour_meta$Name <- paste0(tumour_meta$Pat_ID,"_T")  
normal_meta$Biological <- "Normal"
normal_meta$Name <- paste0(normal_meta$Pat_ID,"_N")  

full_meta <- rbind(normal_meta, tumour_meta) %>% as.data.frame()

columns_to_remove <- c("domain", "phylum", "class", "order", "family", "species")

nejman_data <- nejman_data %>% select(-columns_to_remove)
nejman_data <- nejman_data %>% group_by(genus) %>% dplyr::summarise(across(everything(), sum))
nejman_data <- nejman_data %>% as.data.frame()

normal_data <- nejman_data %>% select(normal_samples, "genus")
normal_data <- normal_data %>% filter(genus %in% total_genus)
normal_otus <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/normal_otu_5percent.csv") %>% as.data.frame()
normal_otus <- normal_otus$x
normal_data <- normal_data %>% filter(genus %in% normal_otu)
rownames(normal_data) <- normal_data$genus
normal_data <- normal_data %>% as.data.frame()

tumour_data <- nejman_data %>% select(tumour_samples, "genus")
tumour_otus <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/tumour_otu_5percent.csv") %>% as.data.frame()
tumour_otus <- tumour_otus$x
tumour_data <- tumour_data %>% filter(genus %in% tumour_otu)
rownames(tumour_data) <- tumour_data$genus
tumour_data <- tumour_data %>% as.data.frame()

full_data <- merge(normal_data, tumour_data, by = "genus", all = TRUE)
full_data[is.na(full_data)] <- 0
rownames(full_data) <- full_data$genus
#genus_to_remove <- c("Diaphorobacter")
full_data <- full_data %>% filter(!genus %in% genus_to_remove)
full_data <- full_data %>% select(!genus)
full_data <- round(full_data)
full_data[full_data < 10] <- 0
full_data <- full_data + 1

full_meta <- full_meta %>% select(Pat_ID, Biological)
##### Step 2: DESeq2 differential abundance #####
dds <- DESeqDataSetFromMatrix(countData = full_data, 
                              colData = full_meta, 
                              design =  ~ Biological)


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
deseq2_results$name <- rownames(full_data)


genus_of_interest <- c("Pseudomonas", "unidentified_Enterobacteriaceae", "Streptococcus", "Tepidiphilus", "Rothia", "Cloacibacterium", 
                       "Klebsiella", "unidentified_Microbacteriaceae", "Aerococcus", "Flaviflexus", "Thauera", "Parapusillimonas", "Citrobacter", 
                       "Diaphorobacter", "Castellaniella")
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


p <- ggplot(data=deseq2_results, aes(x=foldchange, y=-1*log10(pvalue), col=diffexpressed, label = tumour_enriched)) + geom_point(size = 5, alpha = 0.6) + theme_classic() + theme(axis.title.x = element_text(size = rel(2)), 
                                                                                                                                                                                  axis.title.y = element_text(size = rel(2)), 
                                                                                                                                                                                  legend.position="none", 
                                                                                                                                                                                  axis.text.y = element_text(face="bold",size=14),
                                                                                                                                                                                  axis.text.x = element_text(face="bold",size=14)) + scale_color_manual(values = c("deeppink", "darkorange", "deepskyblue3"))
p + geom_label_repel(aes(label=ifelse(diffexpressed == "UP", as.character(name), '')), box.padding= 0.35, point.padding = 0.9,segment.color = 'grey50', max.overlaps = 15, color = "black") + xlab("Fold Change (Tumour vs Normal)") 




EnhancedVolcano(deseq2_results,
                x = 'foldchange',
                y = 'pvalue', 
                lab = deseq2_results$name,
                pCutoff = 1.3, 
                FCcutoff = 0.5, 
                labSize = 2, 
                ylim = c(0, 4))






  
  
