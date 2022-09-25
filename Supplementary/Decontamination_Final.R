# Microbial landscape of NSCLC - Decontamination script using DECONTAM and the negative controls 
# Author - Krupa Thakkar 
# Date - 25th September 2022 

# Loading libraries and packages 
library(tidyverse)
library("phyloseq")
library(ggplot2)
library(decontam)
library(reshape2)

# Loading data 
meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/exp_meta.csv", header = TRUE) %>% as.data.frame()
batch1 <- read.delim("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Raw/batch1_22062021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch1) <- batch1$`#OTU ID`
batch1 <- batch1 %>% select(-`#OTU ID`)
batch2.1 <- read.delim("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Raw/batch2_04112021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch2.1) <- batch2.1$`#OTU ID`
batch2.1 <- batch2.1 %>% select(-`#OTU ID`)
batch2.2 <- read.delim("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Raw/batch2_05112021.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch2.2) <- batch2.2$`#OTU ID`
batch2.2 <- batch2.2 %>% select(-`#OTU ID`)
batch3 <- read.delim("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Raw/batch3_03012022.tsv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(batch3) <- batch3$`#OTU ID`
batch3 <- batch3 %>% select(-`#OTU ID`)

##### Step 1 : Reformatting the meta data and only keeping samples and negative control data + creating a large taxonomy matrix #####
meta_data <- meta_data %>% filter(!Type == "Positive_Control")
meta_data$Sample_or_Control <- NA

for (i in 1:nrow(meta_data)) {
  tmp <- meta_data[i, ] %>% as.data.frame()
  type <- tmp$Type %>% as.character() 
  if (type == "Negative_Control") {
    meta_data[i, "Sample_or_Control"] <- TRUE
  } else {
    meta_data[i, "Sample_or_Control"] <- FALSE
  }
} #TRUE if control and FALSE if sample 
meta_data <- meta_data %>% as.data.frame()

batch1_meta <- meta_data %>% filter(Batch == 1.0)
rownames(batch1_meta) <- batch1_meta$Sample
BATCH1_META <- sample_data(batch1_meta)
batch1_fullnames<- batch1_meta$Sample
batch1_full <- batch1 %>% select(batch1_fullnames)
BATCH1_FULL_OTU <- otu_table(batch1_full, taxa_are_rows = TRUE)

batch2.1_meta <- meta_data %>% filter(Batch == 2.1)
rownames(batch2.1_meta) <- batch2.1_meta$Sample
BATCH2.1_META <- sample_data(batch2.1_meta)
batch2.1_fullnames<- batch2.1_meta$Sample
batch2.1_full <- batch2.1 %>% select(batch2.1_fullnames)
BATCH2.1_FULL_OTU <- otu_table(batch2.1_full, taxa_are_rows = TRUE)

batch2.2_meta <- meta_data %>% filter(Batch == 2.2)
rownames(batch2.2_meta) <- batch2.2_meta$Sample
BATCH2.2_META <- sample_data(batch2.2_meta)
batch2.2_fullnames<- batch2.2_meta$Sample
batch2.2_full <- batch2.2 %>% select(batch2.2_fullnames)
BATCH2.2_FULL_OTU <- otu_table(batch2.2_full, taxa_are_rows = TRUE)

batch3_meta <- meta_data %>% filter(Batch == 3)
rownames(batch3_meta) <- batch3_meta$Sample
BATCH3_META <- sample_data(batch3_meta)
batch3_fullnames<- batch3_meta$Sample
batch3_full <- batch3 %>% select(batch3_fullnames)
BATCH3_FULL_OTU <- otu_table(batch3_full, taxa_are_rows = TRUE)

full_taxa <- c(rownames(batch1), rownames(batch2.1), rownames(batch2.2), rownames(batch3))
full_taxa <- unique(full_taxa)
taxmat <- matrix(nrow = 1789, ncol = 6)
rownames(taxmat) <- full_taxa
colnames(taxmat) <- c("domain", "phylum", "class", "order", "family", "genus")
# Kingdom
for (i in full_taxa) {
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
for (i in full_taxa) {
  full <- i %>% as.character()
  phylum_split <- strsplit(i, split = "p__")
  phylum_split <- phylum_split[[1]][2]
  phylum_split_pt2 <- strsplit(phylum_split, split = ";")
  phylum <- phylum_split_pt2[[1]][1]
  taxmat[i, "phylum"] <- phylum
}
# Class
for (i in full_taxa) {
  full <- i %>% as.character()
  class_split <- strsplit(i, split = "c__")
  class_split <- class_split[[1]][2]
  class_split_pt2 <- strsplit(class_split, split = ";")
  class <- class_split_pt2[[1]][1]
  taxmat[i, "class"] <- class
}
# Order
for (i in full_taxa) {
  full <- i %>% as.character()
  order_split <- strsplit(i, split = "o__")
  order_split <- order_split[[1]][2]
  order_split_pt2 <- strsplit(order_split, split = ";")
  order <- order_split_pt2[[1]][1]
  taxmat[i, "order"] <- order
}
# Family
for (i in full_taxa) {
  full <- i %>% as.character()
  family_split <- strsplit(i, split = "f__")
  family_split <- family_split[[1]][2]
  family_split_pt2 <- strsplit(family_split, split = ";")
  family <- family_split_pt2[[1]][1]
  taxmat[i, "family"] <- family
}
# Genus
for (i in full_taxa) {
  full <- i %>% as.character()
  genus_split <- strsplit(i, split = "g__")
  genus_split <- genus_split[[1]][2]
  genus_split_pt2 <- strsplit(genus_split, split = ";s")
  genus <- genus_split_pt2[[1]][1]
  taxmat[i, "genus"] <- genus
}
TAX <- tax_table(taxmat)

##### Step 2 : Running DECONTAM at a 0.5% prevalence level and saving the outputs #####
phyloseq_object <- phyloseq(BATCH2.2_FULL_OTU, TAX, BATCH2.2_META)
sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$Sample_or_Control == "Negative_Control"
contamdf.prev <- isContaminant(phyloseq_object, method="prevalence", neg="Sample_or_Control", threshold = 0.5)
#saveRDS(contamdf.prev, "~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/batch22_prev_05.rds")

##### Step 3 : Building the final contaminant lists to include Eisenhofer and Slater et al blacklist of contaminants #####
batch1_decontam <- readRDS("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/batch1_prev_05.rds") %>% as.data.frame()
batch21_decontam <- readRDS("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/batch21_prev_05.rds") %>% as.data.frame()
batch22_decontam <- readRDS("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/batch22_prev_05.rds") %>% as.data.frame()
batch3_decontam <- readRDS("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/batch3_prev_05.rds") %>% as.data.frame()
taxmat <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/tax_mat.csv") %>% as.data.frame()
literature_contaminants <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Decontamination/literature_contaminants.csv", header = TRUE) %>% as.data.frame()

rownames(taxmat) <- taxmat$X
taxmat <- taxmat %>% select(-X)
taxmat$name <- NA
all_taxa <- rownames(taxmat)

for (i in rownames(taxmat)) {
  i <- i %>% as.character()
  tmp <- subset(taxmat, rownames(taxmat) == i) %>% as.data.frame()
  genus <- tmp$genus %>% as.character()
  if (is.na(genus) == TRUE) {
    family <- tmp$family %>% as.character()
    taxmat[i, "name"] <- paste0("unidentified_", family)
    if (is.na(family) == TRUE) {
      taxmat[i, "name"] <- NA
    }
  } else {
    taxmat[i, "name"] <- genus
  }
}

taxmat <- taxmat %>% as.data.frame()
uncultured <- "uncultured"

for (i in rownames(taxmat)) {
  i <- i %>% as.character()
  tmp <- subset(taxmat, rownames(taxmat) == i)
  name <- tmp$name %>% as.character()
  if (name %in% uncultured) {
    tmp <- subset(taxmat, rownames(taxmat) == i) %>% as.data.frame()
    family <- tmp$family %>% as.character()
    taxmat[i, "name"] <- paste0("unidentified_", family)
  } 
  else {
    taxmat[i, "name"] <- name
  }
}

literature_contaminants <- literature_contaminants$Contaminant %>% unique()

batch1_decontam <- batch1_decontam %>% filter(contaminant== TRUE)
batch1_decontam$name <- NA
for (i in rownames(batch1_decontam)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch1_decontam[i, "name"] <- name
}
batch1_contaminants <- c(batch1_decontam$name, literature_contaminants)

batch21_decontam <- batch21_decontam %>% filter(contaminant== TRUE)
batch21_decontam$name <- NA
for (i in rownames(batch21_decontam)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch21_decontam[i, "name"] <- name
}
batch21_contaminants <- c(batch21_decontam$name, literature_contaminants)

batch22_decontam <- batch22_decontam %>% filter(contaminant== TRUE)
batch22_decontam$name <- NA
for (i in rownames(batch22_decontam)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch22_decontam[i, "name"] <- name
}
batch22_contaminants <- c(batch22_decontam$name, literature_contaminants)


batch3_decontam <- batch3_decontam %>% filter(contaminant== TRUE)
batch3_decontam$name <- NA
for (i in rownames(batch3_decontam)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch3_decontam[i, "name"] <- name
}
batch3_contaminants <- c(batch3_decontam$name, literature_contaminants)

genus_to_keep <- c("Staphylococcus", "Streptococcus", "Sphingobium", "Roseomonas", "Rhodococcus", "Rothia", "Pseudomonas", "Prevotella", "Neisseria", "Veillonella", "Corynebacterium", "Brevibacterium", "Clostridium", "Escherichia", "Collinsella", "Escherichia-Shigella", "Akkermansia")

batch1_contaminants <- batch1_contaminants %>% unique()
batch1_contaminants <- batch1_contaminants[!(batch1_contaminants %in% genus_to_keep)]

batch21_contaminants <- batch21_contaminants %>% unique()
batch21_contaminants <- batch21_contaminants[!(batch21_contaminants %in% genus_to_keep)]

batch22_contaminants <- batch22_contaminants %>% unique()
batch22_contaminants <- batch22_contaminants[!(batch22_contaminants %in% genus_to_keep)]

batch3_contaminants <- batch3_contaminants %>% unique()
batch3_contaminants <- batch3_contaminants[!(batch3_contaminants %in% genus_to_keep)]

##### Step 4 : Making the final dataframes, calculating % genus removed and saving final output #####
batch1_new <- batch1 %>% as.data.frame()
batch1_new$name <- NA
for (i in rownames(batch1_new)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch1_new[i, "name"] <- name
}
batch1_new <- batch1_new %>% drop_na(name)
batch1_new <- batch1_new %>% group_by(name) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
batch1_new <- batch1_new %>% filter(!name %in% batch1_contaminants)
rownames(batch1_new) <- NULL

batch21_new <- batch2.1 %>% as.data.frame()
batch21_new$name <- NA
for (i in rownames(batch21_new)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch21_new[i, "name"] <- name
}
batch21_new <- batch21_new %>% drop_na(name)
batch21_new <- batch21_new %>% group_by(name) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
batch21_new <- batch21_new %>% filter(!name %in% batch21_contaminants)
rownames(batch21_new) <- NULL

batch22_new <- batch2.2 %>% as.data.frame()
batch22_new$name <- NA
for (i in rownames(batch22_new)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch22_new[i, "name"] <- name
}
batch22_new <- batch22_new %>% drop_na(name)
batch22_new <- batch22_new %>% group_by(name) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
batch22_new <- batch22_new %>% filter(!name %in% batch22_contaminants)
rownames(batch22_new) <- NULL

batch3_new <- batch3 %>% as.data.frame()
batch3_new$name <- NA
for (i in rownames(batch3_new)) {
  i <- i %>% as.character()
  tmp <- taxmat[i, ] %>% as.data.frame()
  name <- tmp$name %>% as.character()
  batch3_new[i, "name"] <- name
}
batch3_new <- batch3_new %>% drop_na(name)
batch3_new <- batch3_new %>% group_by(name) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
batch3_new <- batch3_new %>% filter(!name %in% batch3_contaminants)
rownames(batch3_new) <- NULL

full_data_1 <- merge(batch1_new, batch21_new, by = "name", all = T)
full_data_2 <- merge(full_data_1, batch22_new, by = "name", all = T)
full_data <- merge(full_data_2, batch3_new, by = "name", all = T)

rownames(full_data) <- full_data$name
full_data <- full_data %>% select(-name)
full_data <- full_data %>% as.data.frame()
full_data[is.na(full_data)] <- 0

##### Step 5 : Percentage of spike-in concentrations in the negative controls #####
negative_controls <- meta_data %>% filter(Type == "Negative_Control")
negative_controls$unique_genus <- NA
negative_controls$perc_spikein <- NA
negative_control_samples <- negative_controls$Sample

for (i in 1:nrow(negative_controls)) {
  tmp <- negative_controls[i, ] %>% as.data.frame()
  sample <- tmp$Sample %>% as.character()
  tmp <- full_data %>% select(sample) %>% as.data.frame()
  tmp[tmp < 10] <- NA
  tmp <- tmp %>% na.omit() %>% as.data.frame()
  unique_genus <- nrow(tmp)
  negative_controls[i, "unique_genus"] <- unique_genus
  
}

relative_data <- full_data %>% select(negative_control_samples)
relative_data <- relative_data %>% t() %>% as.data.frame()
relative_data$total <- rowSums(relative_data)
relative_data <- relative_data/relative_data$total
relative_data <- relative_data %>% as.data.frame()
relative_data <- relative_data %>% select(-total)
batch1_neg <- negative_controls %>% filter(Batch == 1)
batch1_neg <- batch1_neg$Sample
batch21_neg <- negative_controls %>% filter(Batch == 2.1)
batch21_neg <- batch21_neg$Sample
batch22_neg <- negative_controls %>% filter(Batch == 2.2)
batch22_neg <- batch22_neg$Sample
batch3_neg <- negative_controls %>% filter(Batch == 3)
batch3_neg <- batch3_neg$Sample

staphy_spikes <- c(batch1_neg, batch21_neg, batch3_neg)

for (i in 1:nrow(negative_controls)) {
  tmp <- negative_controls[i, ] %>% as.data.frame()
  sample <- tmp$Sample %>% as.character()
  
  if(sample %in% staphy_spikes) {
    tmp <- relative_data[sample, ] %>% as.data.frame()
    tmp <- tmp %>% select("Staphylococcus")
    rel_tmp <- tmp$Staphylococcus * 100
  } else {
    tmp <- relative_data[sample, ] %>% as.data.frame()
    tmp <- tmp %>% select("Escherichia-Shigella")
    rel_tmp <- tmp$"Escherichia-Shigella" * 100
  }
  
  negative_controls[i, "perc_spikein"] <- rel_tmp
}

negative_controls <- negative_controls %>% as.data.frame()
negative_controls$Batch <- as.character(negative_controls$Batch)


p <- ggplot(negative_controls, aes(x = Batch, y = perc_spikein)) + geom_boxplot() + theme_classic() + ylab("Relative Abundance of Spike-in in Negative Controls")
p
