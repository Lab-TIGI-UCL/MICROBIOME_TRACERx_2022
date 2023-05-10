# Microbiome in TRACERx NSCLC - Danaher immune scores with tumour enriched bacteria 
# Author - Krupa Thakkar 
# Date- 1st February 2023 

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(vegan)
library(broom)
library(nlme)

# Loading data 
absolute_counts <- read.csv("~/Files/normalised_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Files/Meta_data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

immune_scores <- read.csv("~/Files/Meta_data/immune_Scores.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()

##########
# Step 1 : Reformatting the data 
##########
immune_scores$sample_name <- NA

for (i in 1:nrow(immune_scores)) {
  tmp <- immune_scores[i, ] %>% as.data.frame()
  full_name <- tmp$sample_name_hash %>% as.character()
  split_1 <- strsplit(full_name, split = "--")
  split_name <- split_1[[1]][1]
  new_name <- gsub("X0", "X", split_name) %>% as.character()
  new_name <- gsub("N01", "N", new_name) %>% as.character()
  immune_scores[i, "sample_name"] <- new_name
}

immune_scores <- immune_scores %>% as.data.frame()

meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample

absolute_counts <- absolute_counts %>% select(samples_to_keep)

##########
# Step 2 : Create a larger dataframe with immune scores, diversity and genus of interest 
##########
absolute_counts[absolute_counts < 10] <- 0
absolute_diversity <- absolute_counts %>% t() %>% as.data.frame()
absolute_diversity <- diversity(absolute_diversity, index = "shannon") %>% as.data.frame()
absolute_diversity[absolute_diversity == Inf] <- 0
colnames(absolute_diversity) <- "diversity"
absolute_diversity$danaher_true <- NA

danaher_scores <- immune_scores$sample_name

for (i in rownames(absolute_diversity)) {
  i <- i %>% as.character()
  new_name <- substring(i, 3) %>% as.character()
  
  if (new_name %in% danaher_scores) {
    absolute_diversity[i, "danaher_true"] <- TRUE
    
    immune_tmp <- immune_scores %>% filter(sample_name == new_name) %>% as.data.frame()
    immune_tmp <- immune_tmp %>% select(score, cell_type)
    
    treg <- immune_tmp %>% filter(cell_type == "treg")
    treg <- treg$score
    absolute_diversity[i, "treg"] <- treg 
    
    total_tils <- immune_tmp %>% filter(cell_type == "total_tils")
    total_tils <- total_tils$score
    absolute_diversity[i, "total_tils"] <- total_tils
    
    t_cells <- immune_tmp %>% filter(cell_type == "t_cells")
    t_cells <- t_cells$score 
    absolute_diversity[i, "t_cells"] <- t_cells
    
    th1_cells <- immune_tmp %>% filter(cell_type == "th1_cells")
    th1_cells <- th1_cells$score 
    absolute_diversity[i, "th1_cells"] <- th1_cells
    
    neutrophils <- immune_tmp %>% filter(cell_type == "neutrophils")
    neutrophils <- neutrophils$score 
    absolute_diversity[i, "neutrophils"] <- neutrophils
    
    nk_cells <- immune_tmp %>% filter(cell_type == "nk_cells")
    nk_cells <- nk_cells$score 
    absolute_diversity[i, "nk_cells"] <- nk_cells
    
    macrophages <- immune_tmp %>% filter(cell_type == "macrophages")
    macrophages <- macrophages$score 
    absolute_diversity[i, "macrophages"] <- macrophages
    
    exhausted_cd8_cells <- immune_tmp %>% filter(cell_type == "exhausted_cd8_cells")
    exhausted_cd8_cells <- exhausted_cd8_cells$score
    absolute_diversity[i, "exhausted_cd8_cells"] <- exhausted_cd8_cells
    
    dendritic_cells <- immune_tmp %>% filter(cell_type == "dendritic_cells")
    dendritic_cells <- dendritic_cells$score
    absolute_diversity[i, "dendritic_cells"] <- dendritic_cells
    
    cytotoxic_cells <- immune_tmp %>% filter(cell_type == "cytotoxic_cells")
    cytotoxic_cells <- cytotoxic_cells$score 
    absolute_diversity[i, "cytotoxic_cells"] <- cytotoxic_cells
    
    cd8_cells <- immune_tmp %>% filter(cell_type == "cd8_cells")
    cd8_cells <- cd8_cells$score
    absolute_diversity[i, "cd8_cells"] <- cd8_cells
    
    cd45_cells <- immune_tmp %>% filter(cell_type == "cd45_cells")
    cd45_cells <- cd45_cells$score 
    absolute_diversity[i, "cd45_cells"] <- cd45_cells
    
    cd4_cells <- 	immune_tmp %>% filter(cell_type == "cd4_cells")
    cd4_cells <- cd4_cells$score 
    absolute_diversity[i, "cd4_cells"] <- cd4_cells
    
    b_cells <- immune_tmp %>% filter(cell_type == "b_cells")
    b_cells <- b_cells$score 
    absolute_diversity[i, "b_cells"] <- b_cells
    
    load_tmp <- absolute_counts %>% select(i) %>% as.data.frame()
    load_tmp <- colSums(load_tmp)
    absolute_diversity[i, "microbial_load"] <- load_tmp
    
  }
  else {
    absolute_diversity[i, "danaher_true"] <- FALSE
    
    load_tmp <- absolute_counts %>% select(i) %>% as.data.frame()
    load_tmp <- colSums(load_tmp)
    absolute_diversity[i, "microbial_load"] <- load_tmp
    
  }
}

absolute_diversity <- absolute_diversity %>% as.data.frame()
absolute_diversity$name <- rownames(absolute_diversity)

otu_of_interest <- c("Acidovorax", "Aerococcus", "Brevundimonas", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", "Exiguobacterium", "Flaviflexus", "Fusobacterium", "Hymenobacter", "Lachnoanaerobaculum", "Lactococcus", "Parapusillimonas", "Providencia", "Pseudomonas", "Rothia", "Streptococcus", "Tepidiphilus", "Thauera", "unidentified_Enterobacteriaceae", "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae")
absolute_counts <- subset(absolute_counts, rownames(absolute_counts) %in% otu_of_interest)

absolute_diversity$name <- rownames(absolute_diversity)
absolute_counts <- absolute_counts %>% t() %>% as.data.frame()
absolute_counts$name <- rownames(absolute_counts)

full_dat <- merge(absolute_counts, absolute_diversity, by = "name")
full_dat <- full_dat %>% as.data.frame()

full_dat <- full_dat %>% filter(danaher_true == TRUE) %>% as.data.frame()
full_dat$macrophage_tcell_ratio <- full_dat$macrophage/full_dat$t_cells
full_dat$treg_tcell_ratio <- full_dat$treg/full_dat$t_cells
full_dat <- full_dat %>% as.data.frame()

full_dat$patient_name <- NA

for (i in 1:nrow(full_dat)) {
  
  tmp <- full_dat[i, ] %>% as.data.frame()
  sample <- tmp$name
  tmp_meta <- meta_data %>% filter(Sample == sample)
  pat <- tmp_meta$Patient
  full_dat[i, "patient_name"] <- pat
  
}

##########
# Step 3 : Interaction between immmune cell subsets and microbes
##########
immune_interaction_mat <- matrix(ncol = 14, nrow = 23)
rownames(immune_interaction_mat) <- otu_of_interest
immune_cells <- c("treg", "total_tils", "t_cells", "th1_cells", "neutrophils", "nk_cells", "macrophages", "exhausted_cd8_cells", "cytotoxic_cells", "cd8_cells", "cd45_cells", "cd4_cells", "macrophage_tcell_ratio", "treg_tcell_ratio")
colnames(immune_interaction_mat) <- immune_cells
immune_interaction_mat <- immune_interaction_mat %>% as.data.frame()

for (i in otu_of_interest) {
  i <- i %>% as.character()
  
  for (j in immune_cells) {
    
    if (i == "Providencia") {
      immune_interaction_mat[i, j] <- NA
    } else {
    tmp <- full_dat %>% select(i, j, "patient_name")
    colnames(tmp) <- c("bacteria", "mutation", "patient")
    #tmp$bacteria <- as.factor(as.numeric(tmp$bacteria))
    
    fit <- summary(lme(bacteria ~ mutation, random = ~1|patient, data = tmp))
    # p_val <- fit$tTable[2, 5]
    # 
    # immune_interaction_mat[i, j] <- p_val
    
    tmp <- cov2cor(vcov(fit))
    val <- tmp[1, 2]
    
    immune_interaction_mat[i, j] <- val
    
  }}
}

immune_interaction_mat <- immune_interaction_mat %>% as.data.frame()
immune_interaction_mat <- immune_interaction_mat[complete.cases(immune_interaction_mat), ]
immune_interaction_mat[immune_interaction_mat > 6] <- 6
immune_interaction_mat$genus <- rownames(immune_interaction_mat)
tmp <- melt(immune_interaction_mat)
tmp <- tmp %>% as.data.frame()
tmp$value <- p.adjust(tmp$value, method = "BH")
tmp$value <- -log10(tmp$value)

p <- ggplot(tmp, aes(x = variable, y = genus, fill = value)) + geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 1.30103, name = "-Log10(P Vlaue)") 
p + theme_classic() + xlab("Immune Subsets") + ylab("Genus") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dev.off()
