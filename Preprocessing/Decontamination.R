# Microbiome in NSCLC - Decontamination of 16S data 
# Author - Krupa Thakkar 
# Date - 30th January 2023 

# Loading libraries 
library(tidyverse)
library(decontam)
library(ggplot2)
library(reshape2)
library(phyloseq)

# Loading data + quick formatting 
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

##########
# Step 1 : Data formatting to only keep biological sampels and negative control data + creating a large taxonomy matrix 
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

TAX <- tax_table(taxmat) # Needs to be kept as a matrix and made into a phyloseq object for ease with Decontam 
taxmat <- taxmat %>% as.data.frame() # Changed to a dataframe after being stored as phyloseq object 

meta_data <- meta_data %>% filter(!Type == "Positive_Control") # Removing all positive control samples since only negative controls are of interest 
meta_data$Sample_or_Control <- NA

for (i in 1:nrow(meta_data)) {
  tmp <- meta_data[i, ] %>% as.data.frame()
  type <- tmp$Type %>% as.character() 
  if (type == "Negative_Control") {
    meta_data[i, "Sample_or_Control"] <- "TRUE"
  } else {
    meta_data[i, "Sample_or_Control"] <- "FALSE"
  }
} 
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
batch2.1_full <- batch21 %>% select(batch2.1_fullnames)
BATCH2.1_FULL_OTU <- otu_table(batch2.1_full, taxa_are_rows = TRUE)

batch2.2_meta <- meta_data %>% filter(Batch == 2.2)
rownames(batch2.2_meta) <- batch2.2_meta$Sample
BATCH2.2_META <- sample_data(batch2.2_meta)
batch2.2_fullnames<- batch2.2_meta$Sample
batch2.2_full <- batch22 %>% select(batch2.2_fullnames)
BATCH2.2_FULL_OTU <- otu_table(batch2.2_full, taxa_are_rows = TRUE)

batch3_meta <- meta_data %>% filter(Batch == 3)
rownames(batch3_meta) <- batch3_meta$Sample
BATCH3_META <- sample_data(batch3_meta)
batch3_fullnames<- batch3_meta$Sample
batch3_full <- batch3 %>% select(batch3_fullnames)
BATCH3_FULL_OTU <- otu_table(batch3_full, taxa_are_rows = TRUE)

##########
# Step 2 : Running decontam at a level of 0.5 (most stringent)
##########
phyloseq_object <- phyloseq(BATCH3_FULL_OTU, TAX, BATCH3_META)
sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$Type == "Negative_Control"
contamdf.prev <- isContaminant(phyloseq_object, method="prevalence", neg="is.neg", threshold = 0.5)
write.csv(contamdf.prev, "~/Decontamination/batch3_decontam.csv")


