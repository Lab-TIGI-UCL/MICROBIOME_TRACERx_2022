# Microbial landscape of NSCLC - Tumour enriched bacteria + Ubiquitous bar plot 
# Author - Krupa Thakkar 
# Date - 26th September 2022 

# Loading packges and libraries 
library(tidyverse)
library(ggplot2)
library(reshape2)
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

##### Step 1 : Splitting the meta_data and two dataframes #####
meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(!Sample == "M_LTX258_DNA-BS_GL")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

otu_tables <- rownames(absolute_counts) %>% as.data.frame()
colnames(otu_tables) <- "full_name"
otu_tables$value <- rownames(otu_tables)
otu_tables$value<- sub("^", "otu_", otu_tables$value)

normal_meta <- meta_data %>% filter(Biological == "Normal")
normal_meta <- normal_meta %>% filter(!Patient == "U_LTX745")
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
#absolute_dat_paired[absolute_dat_paired < 10] <- 0
absolute_dat_paired <- absolute_dat_paired + 1

##### Step 2: DESeq2 differential abundance #####
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
                pCutoff = 1.30102999566, 
                FCcutoff = 0.6, 
                labSize = 3,
                pointSize = 5,
                selectLab = c("Aerococcus", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", "Exiguobacterium", "Flaviflexus", 
                              "Hymenobacter", "Klebsiella", "Lachnoanaerobaculum", "Lactococcus", "Parapusillimonas", "Providencia", 
                              "Pseudomonas", "Rothia", "Sphingobium", "Streptococcus", "Tepidiphilus", "Thauera", "unidentified_Enterobacteriaceae", 
                              "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae"),
                labCol = 'black',
                labFace = 'bold', 
                boxedLabels = TRUE, 
                legendPosition = 'bottom', 
                legendLabSize = 10, 
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                ylim = c(0, 3.5), 
                xlim = c(-1.5,1.5), 
                col=c('deeppink', 'deeppink', 'deeppink', 'deepskyblue4'), 
                gridlines.major = FALSE, 
                gridlines.minor = FALSE
)


##### Step 3 : Ubqiuitous boxplots  #####
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()
absolute_counts[absolute_counts < 10] <- 0

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% select(samples_to_keep) %>% as.data.frame()

patient_list <- genus_prevalence$patient %>% unique()

clonality_mat <- matrix(nrow = 22, ncol = 4) %>% as.data.frame()
otu_of_interest <- c("Aerococcus", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", "Exiguobacterium", "Flaviflexus", 
                     "Hymenobacter", "Klebsiella", "Lachnoanaerobaculum", "Lactococcus", "Parapusillimonas", "Providencia", 
                     "Pseudomonas", "Rothia", "Sphingobium", "Streptococcus", "Tepidiphilus", "Thauera","unidentified_Enterobacteriaceae", 
                     "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae")
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
clonality_mat <- matrix(nrow = 105, ncol = 4) %>% as.data.frame()
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
