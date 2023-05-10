# Microbial landscape of NSCLC in TRACERx - Full preprocessing (Data lock)
# Author - Krupa Thakkar 
# Date - 10th March 2023 

# Loading libraries 
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(nlme)
library('dplyr')
library(edgeR)
library(plyr)
require(reshape2)
library(gridExtra)
library(scales)
library("readr")
library(limma)
library(snm)
library(vegan)

# Data lock definition script 
meta_data <- read.csv("~/Meta_data/exp_meta.csv", header = TRUE) %>% as.data.frame()

batch1 <- read.delim("~/Raw/batch1_22062021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch1) <- batch1$`#OTU ID`
batch1 <- batch1 %>% select(-`#OTU ID`)

batch21 <- read.delim("~/Raw/batch2_04112021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch21) <- batch21$`#OTU ID`
batch21 <- batch21 %>% select(-`#OTU ID`)

batch22 <- read.delim("~/Raw/batch2_05112021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch22) <- batch22$`#OTU ID`
batch22 <- batch22 %>% select(-`#OTU ID`)

batch3 <- read.delim("~/Raw/batch3_03012022.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch3) <- batch3$`#OTU ID`
batch3 <- batch3 %>% select(-`#OTU ID`)

batch1_contaminants <- read.csv("~/Decontamination/batch1_decontam.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch1_contaminants) <- batch1_contaminants[,1]
batch1_contaminants[,1] <- NULL

batch21_contaminants <- read.csv("~/Decontamination/batch21_decontam.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch21_contaminants) <- batch21_contaminants[,1]
batch21_contaminants[,1] <- NULL

batch22_contaminants <- read.csv("~/Decontamination/batch22_decontam.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch22_contaminants) <- batch22_contaminants[,1]
batch22_contaminants[,1] <- NULL

batch3_contaminants <- read.csv("~/Decontamination/batch3_decontam.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch3_contaminants) <- batch3_contaminants[,1]
batch3_contaminants[,1] <- NULL

literature_contaminants <- read_xlsx("~/Files/16S_data/Decontamination/literature_contaminants_new.xlsx") %>% as.data.frame()
##########
# Step 1 : Data formatting to only keep biological sampels and negative control data + creating a large taxonomy matrix + keeping everything at a genus level 
##########
all_taxa <- c(rownames(batch1), rownames(batch21), rownames(batch22), rownames(batch3)) %>% unique()
taxmat <- matrix(nrow = 1789, ncol = 6)
rownames(taxmat) <- all_taxa
colnames(taxmat) <- c("domain", "phylum", "class", "order", "family", "genus")
# Kingdom
for (i in all_taxa) {
  full <- i %>% as.character()
  kingdom_split <- strsplit(i, split = "d__")
  kingdom_split <- kingdom_split[[1]][2]
  kingdom_split_pt2 <- strsplit(kingdom_split, split = ";")
  kingdom <- kingdom_split_pt2[[1]][1]
  if (is.na(kingdom) == TRUE) {
    taxmat[i, "domain"] <- NA
  } else {
    taxmat[i, "domain"] <- kingdom
  }}
# Phylum
for (i in all_taxa) {
  full <- i %>% as.character()
  phylum_split <- strsplit(i, split = "p__")
  phylum_split <- phylum_split[[1]][2]
  phylum_split_pt2 <- strsplit(phylum_split, split = ";")
  phylum <- phylum_split_pt2[[1]][1]
  taxmat[i, "phylum"] <- phylum
}
# Class
for (i in all_taxa) {
  full <- i %>% as.character()
  class_split <- strsplit(i, split = "c__")
  class_split <- class_split[[1]][2]
  class_split_pt2 <- strsplit(class_split, split = ";")
  class <- class_split_pt2[[1]][1]
  taxmat[i, "class"] <- class
}
# Order
for (i in all_taxa) {
  full <- i %>% as.character()
  order_split <- strsplit(i, split = "o__")
  order_split <- order_split[[1]][2]
  order_split_pt2 <- strsplit(order_split, split = ";")
  order <- order_split_pt2[[1]][1]
  taxmat[i, "order"] <- order
}
# Family
for (i in all_taxa) {
  full <- i %>% as.character()
  family_split <- strsplit(i, split = "f__")
  family_split <- family_split[[1]][2]
  family_split_pt2 <- strsplit(family_split, split = ";")
  family <- family_split_pt2[[1]][1]
  taxmat[i, "family"] <- family
}
# Genus
for (i in all_taxa) {
  full <- i %>% as.character()
  genus_split <- strsplit(i, split = "g__")
  genus_split <- genus_split[[1]][2]
  genus_split_pt2 <- strsplit(genus_split, split = ";s")
  genus <- genus_split_pt2[[1]][1]
  taxmat[i, "genus"] <- genus
}

taxmat <- taxmat %>% as.data.frame()
taxmat$final_name <- NA # Adding in the new final name which is either the genus or unidentified family 

for (i in rownames(taxmat)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  genus <- tmp_dat$genus %>% as.character()
  
  if (is.na(genus) == TRUE) {
    family <- tmp_dat$family %>% as.character()
    new_name <- paste0("unidentified_",family) %>% as.character()
    taxmat[i, "final_name"] <- new_name
  } else {
    
    taxmat[i, "final_name"] <- genus
  }}

taxmat <- taxmat %>% as.data.frame()

##########
# Step 2 : Building the larger decontamination list 
##########

to_keep_literature <- c("Lactobacillus", "Staphylococcus", "Streptococcus", "Rhodococcus", "Rothia", "Pseudomonas", "Prevotella","Brevundimonas",  "Neisseria", "Veillonella", "Escherichia", "Escherichia-Shigella", "Acidovorax", "Fusobacterium")

all_literature_contaminants <- literature_contaminants$Contaminant %>% unique()

batch1_contam_list <- batch1_contaminants
batch1_contam_list <- batch1_contam_list %>% filter(contaminant == "TRUE") %>% as.data.frame()
batch1_contam_list$final_name <- NA
for (i in rownames(batch1_contam_list)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch1_contam_list[i, "final_name"] <- new_name 
}
batch1_contam_list <- batch1_contam_list %>% as.data.frame()
batch1_contam_list <- batch1_contam_list$final_name 
batch1_contam_list <- c(batch1_contam_list, all_literature_contaminants) %>% unique()
batch1_contam_list <- batch1_contam_list[!(batch1_contam_list %in% to_keep_literature)]


batch21_contam_list <- batch21_contaminants
batch21_contam_list <- batch21_contam_list %>% filter(contaminant == "TRUE") %>% as.data.frame()
batch21_contam_list$final_name <- NA
for (i in rownames(batch21_contam_list)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch21_contam_list[i, "final_name"] <- new_name 
}
batch21_contam_list <- batch21_contam_list %>% as.data.frame()
batch21_contam_list <- batch21_contam_list$final_name 
batch21_contam_list <- c(batch21_contam_list, all_literature_contaminants) %>% unique()
batch21_contam_list <- batch21_contam_list[!(batch21_contam_list %in% to_keep_literature)]

batch22_contam_list <- batch22_contaminants
batch22_contam_list <- batch22_contam_list %>% filter(contaminant == "TRUE") %>% as.data.frame()
batch22_contam_list$final_name <- NA
for (i in rownames(batch22_contam_list)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch22_contam_list[i, "final_name"] <- new_name 
}
batch22_contam_list <- batch22_contam_list %>% as.data.frame()
batch22_contam_list <- batch22_contam_list$final_name 
batch22_contam_list <- c(batch22_contam_list, all_literature_contaminants) %>% unique()
batch22_contam_list <- batch22_contam_list[!(batch22_contam_list %in% to_keep_literature)]


batch3_contam_list <- batch3_contaminants
batch3_contam_list <- batch3_contam_list %>% filter(contaminant == "TRUE") %>% as.data.frame()
batch3_contam_list$final_name <- NA
for (i in rownames(batch3_contam_list)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch3_contam_list[i, "final_name"] <- new_name 
}
batch3_contam_list <- batch3_contam_list %>% as.data.frame()
batch3_contam_list <- batch3_contam_list$final_name 
batch3_contam_list <- c(batch3_contam_list, all_literature_contaminants) %>% unique()
batch3_contam_list <- batch3_contam_list[!(batch3_contam_list %in% to_keep_literature)]

##########
# Step 3 : Removing the contaminants from the list 
##########
batch1 <- batch1 %>% as.data.frame()
batch1$full_name <- NA
for (i in rownames(batch1)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch1[i, "full_name"] <- new_name
}
batch1 <- batch1 %>% as.data.frame()
batch1 <- batch1 %>% filter(!full_name %in% batch1_contam_list)


batch21 <- batch21 %>% as.data.frame()
batch21$full_name <- NA
for (i in rownames(batch21)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch21[i, "full_name"] <- new_name
}
batch21 <- batch21 %>% as.data.frame()
batch21 <- batch21 %>% filter(!full_name %in% batch21_contam_list)

batch22 <- batch22 %>% as.data.frame()
batch22$full_name <- NA
for (i in rownames(batch22)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch22[i, "full_name"] <- new_name
}
batch22 <- batch22 %>% as.data.frame()
batch22 <- batch22 %>% filter(!full_name %in% batch22_contam_list)

batch3 <- batch3 %>% as.data.frame()
batch3$full_name <- NA
for (i in rownames(batch3)) {
  full_name <- i %>% as.character()
  tmp_dat <- taxmat[i, ] %>% as.data.frame()
  new_name <- tmp_dat$final_name %>% as.character()
  batch3[i, "full_name"] <- new_name
}
batch3 <- batch3 %>% as.data.frame()
batch3 <- batch3 %>% filter(!full_name %in% batch3_contam_list)

full_data_1 <- merge(batch1, batch21, by = "full_name", all = T)
full_data_2 <- merge(full_data_1, batch22, by = "full_name", all = T)
full_data <- merge(full_data_2, batch3, by = "full_name", all = T)

rownames(full_data) <- full_data$full_name
full_data <- full_data %>% select(-full_name)
full_data <- full_data %>% as.data.frame()
full_data[is.na(full_data)] <- 0

full_data <- full_data %>% as.data.frame()

samples_to_keep <- meta_data %>% filter(Type == "Sample")
samples_to_keep <- samples_to_keep %>% filter(!Batch == "2.2")
samples_to_keep <- samples_to_keep$Sample

new_data <- full_data %>% select(samples_to_keep)
new_data <- new_data %>% as.data.frame()
new_data$name <- rownames(new_data)
to_remove <- c("Chloroplast", "Mitochondria")
new_data <- new_data %>% filter(!name %in% to_remove)
new_data <- new_data %>% as.data.frame() %>% select(-name)
new_data <- new_data %>% as.data.frame()

full_data <- new_data
meta_data <- read.csv("~/Meta_data/exp_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()

samples_meta <- meta_data %>% filter(Type == "Sample")
samples_meta <- samples_meta %>% filter(!Batch == "2.2")

counts <- full_data
dge <- DGEList(counts = counts)
v <- voom(counts, plot=TRUE, normalize="quantile")
v_matrix <- v$E

rownames(samples_meta) <- samples_meta$Sample
bio.var <- model.matrix(~Biological, data = samples_meta)
adj.var <- model.matrix(~Batch, data = samples_meta)

snm_object <- snm(raw.dat = v_matrix, 
                  bio.var = bio.var, 
                  adj.var = adj.var, 
                  rm.adj = TRUE, 
                  verbose = TRUE, 
                  diagnose = TRUE)
snm_data <- (snm_object$norm.dat)

full_data <- snm_data
