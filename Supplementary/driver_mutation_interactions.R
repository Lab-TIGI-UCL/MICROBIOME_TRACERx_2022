# Microbiome in TRACERx NSCLC - Mutational genotypes with tumour enriched bacteria 
# Author - Krupa Thakkar 
# Date - 2nd February 2023 

# Loading libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
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

mut_table <- get(load(file = "~//Files/421_mut_table/patientMutTable.RData"))
small_mut_table <- mut_table %>% select(SampleID, Hugo_Symbol, DriverMut)
tumour_ID <- get(load(file = "~/Files/421_mut_table/tumour_sample_id_df.RData")) %>% as.data.frame()

##########
# Step 1 : Data formatting 
##########
absolute_counts[absolute_counts < 10] <- 0
#absolute_counts[absolute_counts > 9] <- 1

meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)
otu_of_interest <- c("Acidovorax", "Aerococcus", "Brevundimonas", "Castellaniella", "Citrobacter", "Cloacibacterium", "Diaphorobacter", "Exiguobacterium", "Flaviflexus", "Fusobacterium", "Hymenobacter", "Lachnoanaerobaculum", "Lactococcus", "Parapusillimonas", "Pseudomonas", "Rothia", "Streptococcus", "Tepidiphilus", "Thauera", "unidentified_Enterobacteriaceae", "unidentified_Microbacteriaceae", "unidentified_Micrococcaceae")

##########
# Step 2 : Creating large dataframes for LUAD and LUSC 
##########
luad_mat <- absolute_counts
luad_pat <- meta_data %>% filter(Histology == "LUAD")
luad_mat <- luad_mat %>% select(luad_pat$Sample)
luad_genus <- otu_of_interest
luad_mat <- luad_mat[luad_genus, ]
luad_mat <- luad_mat %>% t() %>% as.data.frame()
luad_mat$KRAS <- NA
luad_mat$EGFR <- NA
luad_mat$TP53 <- NA
luad_mat$KEAP1 <- NA
luad_mat$FAT1 <- NA
luad_mat$patient_name <- NA

for (i in rownames(luad_mat)) {
  sample_name <- meta_data %>% filter(Sample == i)
  sample_name <- sample_name[1, ] %>% as.data.frame()
  sample_name <- sample_name$Sample %>% as.character()
  print(paste0(i, "_sample_name"))
  patient_name <- tumour_ID %>% filter(sample_id == sample_name)
  patient_name <- patient_name[1, ] %>% as.data.frame()
  patient_name <- patient_name$patient_id %>% as.character()
  print(paste0(i, "_patient_name"))
  mut_data <- small_mut_table %>% filter(SampleID == patient_name) %>% as.data.frame()
  
  luad_mat[i, "patient_name"] <- patient_name
  
  KRAS_dat <- mut_data %>% filter(Hugo_Symbol == "KRAS")
  KRAS_dat <- KRAS_dat[1, ] %>% as.data.frame()
  KRAS_dat <- KRAS_dat$DriverMut %>% as.character()
  luad_mat[i, "KRAS"] <- KRAS_dat
  print(paste0(i, "_KRAS_done"))
  
  EGFR_dat <- mut_data %>% filter(Hugo_Symbol == "EGFR")
  EGFR_dat <- EGFR_dat[1, ] %>% as.data.frame()
  EGFR_dat <- EGFR_dat$DriverMut %>% as.character()
  luad_mat[i, "EGFR"] <- EGFR_dat
  print(paste0(i, "_EGFR_done"))
  
  TP53_dat <- mut_data %>% filter(Hugo_Symbol == "TP53")
  TP53_dat <- TP53_dat[1, ] %>% as.data.frame()
  TP53_dat <- TP53_dat$DriverMut %>% as.character()
  luad_mat[i, "TP53"] <- TP53_dat
  print(paste0(i, "_TP53_done"))
  
  FAT1_dat <- mut_data %>% filter(Hugo_Symbol == "FAT1")
  FAT1_dat <- FAT1_dat[1, ] %>% as.data.frame()
  FAT1_dat <- FAT1_dat$DriverMut %>% as.character()
  luad_mat[i, "FAT1"] <- FAT1_dat
  print(paste0(i, "_FAT1_done"))
  
  KEAP1_dat <- mut_data %>% filter(Hugo_Symbol == "KEAP1")
  KEAP1_dat <- KEAP1_dat[1, ] %>% as.data.frame()
  KEAP1_dat <- KEAP1_dat$DriverMut %>% as.character()
  luad_mat[i, "KEAP1"] <- KEAP1_dat
  print(paste0(i, "_KEAP1_done"))
}

lusc_mat <- absolute_counts
lusc_pat <- meta_data %>% filter(Histology == "LUSC")
lusc_mat <- lusc_mat %>% select(lusc_pat$Sample)
lusc_genus <- otu_of_interest
lusc_mat <- lusc_mat[lusc_genus, ]
lusc_mat <- lusc_mat %>% t() %>% as.data.frame()
lusc_mat$TP53 <- NA
lusc_mat$NOTCH1 <- NA
lusc_mat$CDKN2A <- NA
lusc_mat$FAT1 <- NA
lusc_mat$PIK3CA <- NA
lusc_mat$patient_name <- NA

for (i in rownames(lusc_mat)) {
  sample_name <- meta_data %>% filter(Sample == i)
  sample_name <- sample_name[1, ] %>% as.data.frame()
  sample_name <- sample_name$Sample %>% as.character()
  print(paste0(i, "_sample_name"))
  patient_name <- tumour_ID %>% filter(sample_id == sample_name)
  patient_name <- patient_name[1, ] %>% as.data.frame()
  patient_name <- patient_name$patient_id %>% as.character()
  print(paste0(i, "_patient_name"))
  mut_data <- small_mut_table %>% filter(SampleID == patient_name) %>% as.data.frame()
  
  lusc_mat[i, "patient_name"] <- patient_name
  
  TP53_dat <- mut_data %>% filter(Hugo_Symbol == "TP53")
  TP53_dat <- TP53_dat[1, ] %>% as.data.frame()
  TP53_dat <- TP53_dat$DriverMut %>% as.character()
  lusc_mat[i, "TP53"] <- TP53_dat
  print(paste0(i, "_KRAS_done"))
  
  NOTCH1_dat <- mut_data %>% filter(Hugo_Symbol == "NOTCH1")
  NOTCH1_dat <- NOTCH1_dat[1, ] %>% as.data.frame()
  NOTCH1_dat <- NOTCH1_dat$DriverMut %>% as.character()
  lusc_mat[i, "NOTCH1"] <- NOTCH1_dat
  print(paste0(i, "_NOTCH1_done"))
  
  CDKN2A_dat <- mut_data %>% filter(Hugo_Symbol == "CDKN2A")
  CDKN2A_dat <- CDKN2A_dat[1, ] %>% as.data.frame()
  CDKN2A_dat <- CDKN2A_dat$DriverMut %>% as.character()
  lusc_mat[i, "CDKN2A"] <- CDKN2A_dat
  print(paste0(i, "_CDKN2A_done"))
  
  FAT1_dat <- mut_data %>% filter(Hugo_Symbol == "FAT1")
  FAT1_dat <- FAT1_dat[1, ] %>% as.data.frame()
  FAT1_dat <- FAT1_dat$DriverMut %>% as.character()
  lusc_mat[i, "FAT1"] <- FAT1_dat
  print(paste0(i, "_FAT1_done"))
  
  PIK3CA_dat <- mut_data %>% filter(Hugo_Symbol == "TP53")
  PIK3CA_dat <- PIK3CA_dat[1, ] %>% as.data.frame()
  PIK3CA_dat <- PIK3CA_dat$DriverMut %>% as.character()
  lusc_mat[i, "PIK3CA"] <- PIK3CA_dat
  print(paste0(i, "_PIK3CA_done"))
}

##########
# Step 3 : Linear regression analysis 
##########
luad_mat <- luad_mat[!is.na(luad_mat$patient_name),]
luad_mat[is.na(luad_mat)] <- "FALSE"
luad_mat <- luad_mat %>% as.data.frame()
luad_heatmap_mat <- matrix(nrow = 22, ncol = 5)
colnames(luad_heatmap_mat) <- c("TP53", "KRAS", "EGFR", "FAT1", "KEAP1") # For LUAD 
rownames(luad_heatmap_mat) <- otu_of_interest
luad_mutations <- c("TP53", "KRAS", "EGFR", "FAT1", "KEAP1")

for (i in otu_of_interest) {
  i <- i %>% as.character()
  print(i)
  
  for (j in luad_mutations) {
    tmp <- luad_mat %>% select(i, j, "patient_name")
    colnames(tmp) <- c("bacteria", "mutation", "patient")
    #tmp$bacteria <- as.factor(as.numeric(tmp$bacteria))
    vals <- c("FALSE", "TRUE")
    tmp$mutation <- factor(tmp$mutation, levels = vals)
    
    fit <- summary(lme(bacteria ~ mutation, random = ~1|patient, data = tmp))
    p_val <- fit$tTable[2, 5]
    
    luad_heatmap_mat[i, j] <- p_val
    
  }
}

luad_heatmap_mat <- luad_heatmap_mat %>% as.data.frame()
luad_heatmap_mat[luad_heatmap_mat > 6] <- 6
luad_heatmap_mat$genus <- rownames(luad_heatmap_mat)
tmp <- melt(luad_heatmap_mat)
tmp <- tmp %>% as.data.frame()
tmp$value <- p.adjust(tmp$value, method = "BH")
tmp$value <- -log10(tmp$value)

p <- ggplot(tmp, aes(x = variable, y = genus, fill = value)) + geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 1.30103, name = "-Log10(P Vlaue)") 
p + theme_classic() + xlab("Genes") + ylab("Genus")

lusc_mat <- lusc_mat[!is.na(lusc_mat$patient_name),]
lusc_mat[is.na(lusc_mat)] <- "FALSE"
lusc_mat <- lusc_mat %>% as.data.frame()
lusc_heatmap_mat <- matrix(nrow = 22, ncol = 5)
colnames(lusc_heatmap_mat) <- c("TP53", "NOTCH1", "CDKN2A", "FAT1", "PIK3CA") # For LUSC
rownames(lusc_heatmap_mat) <- otu_of_interest
lusc_mutations <- c("TP53", "NOTCH1", "CDKN2A", "FAT1", "PIK3CA")

for (i in otu_of_interest) {
  i <- i %>% as.character()
  print(i)
  
  for (j in lusc_mutations) {
    tmp <- lusc_mat %>% select(i, j, "patient_name")
    colnames(tmp) <- c("bacteria", "mutation", "patient")
    #tmp$bacteria <- as.factor(as.numeric(tmp$bacteria))
    vals <- c("FALSE", "TRUE")
    tmp$mutation <- factor(tmp$mutation, levels = vals)
    
    fit <- summary(lme(bacteria ~ mutation, random = ~1|patient, data = tmp))
    p_val <- fit$tTable[2, 5]
    
    lusc_heatmap_mat[i, j] <- p_val
    
  }
}

lusc_heatmap_mat <- lusc_heatmap_mat %>% as.data.frame()
lusc_heatmap_mat[lusc_heatmap_mat > 6] <- 6
lusc_heatmap_mat$genus <- rownames(lusc_heatmap_mat)
tmp <- melt(lusc_heatmap_mat)
tmp <- tmp %>% as.data.frame()
tmp$value <- p.adjust(tmp$value, method = "BH")
tmp$value <- -log10(tmp$value)

p <- ggplot(tmp, aes(x = variable, y = genus, fill = value)) + geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 1.30103, name = "-Log10(P Vlaue)") 
p + theme_classic() + xlab("Genes") + ylab("Genus")


##########
# Step 4 : Patient specific dataframe instead of region specific with a fishers exact test 
##########
# LUAD
luad_patients <- unique(luad_mat$patient_name)
luad_pat_df <- matrix(nrow = 102, ncol = 27) %>% as.data.frame()
rownames(luad_pat_df) <- luad_patients
colnames(luad_pat_df) <- c(otu_of_interest, luad_mutations)

for (i in luad_patients) {
  i <- i %>% as.character()
  mat <- luad_mat %>% filter(patient_name == i)
  
  for (j in colnames(luad_pat_df)) {
    j <- j %>% as.character()
    
    if (j %in% otu_of_interest) {
      j <- j %>% as.character()
      tmp <- mat %>% select(j) %>% as.data.frame() 
      tmp <- tmp %>% t() %>% as.data.frame()
      tmp$total <- rowSums(tmp)
      if(tmp$total > 0) {
        luad_pat_df[i, j] <- "PRESENT"
      } else {
      luad_pat_df[i, j]  <- "ABSENT"
      
      }
    } else {
      tmp <- mat[1, ]
      tmp <- tmp %>% select(j) %>% as.data.frame()
      colnames(tmp) <- "val"
      
      luad_pat_df[i, j] <- tmp$val 
    }
  }
}

luad_heatmap_mat <- matrix(nrow = 22, ncol = 5)
colnames(luad_heatmap_mat) <- c("TP53", "KRAS", "EGFR", "FAT1", "KEAP1") # For LUAD 
rownames(luad_heatmap_mat) <- otu_of_interest

for(i in rownames(luad_heatmap_mat)) {
  dat_1 <- subset(luad_heatmap_mat, rownames(luad_heatmap_mat) == i)
  mutations <- colnames(dat_1)
  for (j in mutations) {
    full_dat <- luad_pat_df %>% select(i, j) %>% as.data.frame() #change this to LUAD or LUSC
    colnames(full_dat) <- c("Bacteria", "Mutation")
    cont_table <- table(full_dat$Bacteria, full_dat$Mutation)
    
    if (nrow(cont_table) == 1) {
      luad_heatmap_mat[i, j] <- NA
    } else {
    p_val <- fisher.test(cont_table)
    p_val <- p_val$p.value
    # p_val <- p.adjust(p_val)
    # p_val <- -log10(p_val)
    luad_heatmap_mat[i, j] <- p_val}
}}

luad_heatmap_mat <- luad_heatmap_mat %>% as.data.frame()
luad_heatmap_mat[luad_heatmap_mat > 6] <- 6
luad_heatmap_mat$genus <- rownames(luad_heatmap_mat)
tmp <- melt(luad_heatmap_mat)
tmp <- tmp %>% as.data.frame()
tmp$value <- p.adjust(tmp$value, method = "BH")
tmp$value <- -log10(tmp$value)

p <- ggplot(tmp, aes(x = variable, y = genus, fill = value)) + geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 1.30103, name = "-Log10(P Vlaue)") 
p + theme_classic() + xlab("Genes") + ylab("Genus")

# LUSC
lusc_mutations <- c("TP53", "NOTCH1", "CDKN2A", "FAT1", "PIK3CA")
lusc_mat <- lusc_mat[!is.na(lusc_mat$patient_name),]
lusc_mat[is.na(lusc_mat)] <- "FALSE"
lusc_mat <- lusc_mat %>% as.data.frame()
lusc_patients <- unique(lusc_mat$patient_name)
lusc_pat_df <- matrix(nrow = 66, ncol = 27) %>% as.data.frame()
rownames(lusc_pat_df) <- lusc_patients
colnames(lusc_pat_df) <- c(otu_of_interest, lusc_mutations)

for (i in lusc_patients) {
  i <- i %>% as.character()
  mat <- lusc_mat %>% filter(patient_name == i)
  
  for (j in colnames(lusc_pat_df)) {
    j <- j %>% as.character()
    
    if (j %in% otu_of_interest) {
      j <- j %>% as.character()
      tmp <- mat %>% select(j) %>% as.data.frame() 
      tmp <- tmp %>% t() %>% as.data.frame()
      tmp$total <- rowSums(tmp)
      if(tmp$total > 0) {
        lusc_pat_df[i, j] <- "PRESENT"
      } else {
        lusc_pat_df[i, j]  <- "ABSENT"
        
      }
    } else {
      tmp <- mat[1, ]
      tmp <- tmp %>% select(j) %>% as.data.frame()
      colnames(tmp) <- "val"
      
      lusc_pat_df[i, j] <- tmp$val 
    }
  }
}

lusc_heatmap_mat <- matrix(nrow = 22, ncol = 5)
colnames(lusc_heatmap_mat) <- c("TP53", "NOTCH1", "CDKN2A", "FAT1", "PIK3CA") # For lusc 
rownames(lusc_heatmap_mat) <- otu_of_interest

for(i in rownames(lusc_heatmap_mat)) {
  dat_1 <- subset(lusc_heatmap_mat, rownames(lusc_heatmap_mat) == i)
  mutations <- colnames(dat_1)
  for (j in mutations) {
    full_dat <- lusc_pat_df %>% select(i, j) %>% as.data.frame() #change this to lusc or LUSC
    colnames(full_dat) <- c("Bacteria", "Mutation")
    cont_table <- table(full_dat$Bacteria, full_dat$Mutation)
    
    if (nrow(cont_table) == 1) {
      lusc_heatmap_mat[i, j] <- NA
    } else {
      p_val <- fisher.test(cont_table)
      p_val <- p_val$p.value
      # p_val <- p.adjust(p_val)
      # p_val <- -log10(p_val)
      lusc_heatmap_mat[i, j] <- p_val}
  }}

lusc_heatmap_mat <- lusc_heatmap_mat %>% as.data.frame()
lusc_heatmap_mat[lusc_heatmap_mat > 6] <- 6
lusc_heatmap_mat$genus <- rownames(lusc_heatmap_mat)
tmp <- melt(lusc_heatmap_mat)
tmp <- tmp %>% as.data.frame()
tmp$value <- p.adjust(tmp$value, method = "BH")
tmp$value <- -log10(tmp$value)

p <- ggplot(tmp, aes(x = variable, y = genus, fill = value)) + geom_tile() + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 1.30103, name = "-Log10(P Vlaue)") 
p + theme_classic() + xlab("Genes") + ylab("Genus")



