# Microbial landscape of NSCLC - Identifing bacteria associated with lung disorders to validate pipeline 
# Author - Krupa Thakkar  
# Date - 28th February 2023 

# Loading packages 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(vegan)
library(broom)
library(nlme)

# Loading data 
matrix <- read.csv("~/Files/normalised_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(matrix) <- matrix[, 1] 
matrix[, 1] <- NULL
matrix <- matrix %>% as.data.frame()
counts <- matrix

meta_data <- read.csv("~/Files/Meta_data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

##### Section 1 : Reformatting the data #####
batches <- c("2.2")
meta_data <- meta_data%>% filter(!Batch %in% batches)
samples_to_keep <- c("Tumour")
meta_data <- meta_data %>% filter(Biological %in% samples_to_keep)
samples_to_keep <- meta_data$Sample

counts <- counts %>% dplyr::select(samples_to_keep)
counts[counts < 10] <- 0

##### Section 2 : Building a large dataframe #####
full_dat <- colnames(counts) %>% as.data.frame()
colnames(full_dat) <- "samples"
rownames(full_dat) <- full_dat$samples
full_dat$unique_genus <- NA
full_dat$microbial_load <- NA
full_dat$staphylococcus <- NA
full_dat$pseudomonas <- NA
full_dat$streptococcus <- NA
full_dat$total_score <- NA
full_dat$emphysema <- NA
full_dat$copd <- NA

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp <- counts %>% select(i) %>% as.data.frame()
  
  unique_genus <- tmp %>% as.data.frame()
  unique_genus[unique_genus == 0] <- NA
  unique_genus <- unique_genus %>% na.omit()
  unique_genus <- nrow(unique_genus) 
  full_dat[i, "unique_genus"] <- unique_genus
  
  microbial_load <- tmp %>% as.data.frame()
  microbial_load[microbial_load == 0] <- NA
  microbial_load <- microbial_load %>% na.omit()
  microbial_load <- colSums(microbial_load) 
  full_dat[i, "microbial_load"] <- microbial_load
  
  bacteria_tmp <- tmp %>% t() %>% as.data.frame()
  bacteria_tmp <- bacteria_tmp %>% select("Staphylococcus", "Streptococcus", "Pseudomonas")
  staphy <- bacteria_tmp$Staphylococcus
  full_dat[i, "staphylococcus"] <- staphy
  strepto <- bacteria_tmp$Streptococcus
  full_dat[i, "streptococcus"] <- strepto
  pseudo <- bacteria_tmp$Pseudomonas
  full_dat[i, "pseudomonas"] <- pseudo
  
  bacteria_tmp[bacteria_tmp < 10] <- NA
  bacteria_tmp <- bacteria_tmp %>% t() %>% as.data.frame()
  bacteria_tmp <- bacteria_tmp %>% na.omit() %>% as.data.frame()
  bacteria_tmp <- nrow(bacteria_tmp) 
  full_dat[i, "total_score"] <- bacteria_tmp
  
}
full_dat <- full_dat %>% as.data.frame()
full_dat$pathology <- NA

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp_meta <- meta_data %>% filter(Sample == i)
  pathology <- tmp_meta$PathTissueFind 
  full_dat[i, "pathology"] <- pathology
}

# Manually checked all pathology 
#write.csv(full_dat, "~/Files/association_meta_new.csv")

##### Section 3 : Testing associations  #####
full_dat <- read.csv("~/Files/association_meta_new.csv", header = TRUE) %>% as.data.frame() # <- an older one but still can overwrite the load and genus 
full_dat[, 1] <- NULL
rownames(full_dat) <- full_dat[, 1] 
full_dat <- full_dat %>% as.data.frame()
full_dat[is.na(full_dat)] <- FALSE
order <- c("TRUE", "FALSE")
full_dat$pneumonia <- factor(full_dat$pneumonia, 
                             level = order)

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp <- counts %>% select(i) %>% as.data.frame()
  unique_genus <- tmp %>% as.data.frame()
  unique_genus[unique_genus == 0] <- NA
  unique_genus <- unique_genus %>% na.omit()
  unique_genus <- nrow(unique_genus) 
  full_dat[i, "unique_genus"] <- unique_genus
  
  microbial_load <- tmp %>% as.data.frame()
  microbial_load[microbial_load == 0] <- NA
  microbial_load <- microbial_load %>% na.omit()
  microbial_load <- colSums(microbial_load) 
  full_dat[i, "microbial_load"] <- microbial_load
  
  
}

full_dat$patient <- NA

for (i in samples_to_keep) {
  patient <- meta_data %>% filter(Sample == i)
  patient <- patient$Patient %>% as.character()
  
  full_dat[i, "patient"] <- patient
}

fit <- summary(lme(unique_genus ~ pneumonia, random = ~1|patient, data = full_dat))
p_val <- fit$tTable[2, 5]
p_val <- p.adjust(p_val)


p <- ggplot(full_dat, aes(x = pneumonia, y = microbial_load, fill = pneumonia)) + geom_boxplot()
p +  theme_classic() + stat_compare_means() + scale_fill_brewer(palette = "BuPu") + xlab("Pneumonia Diagnosis") + ylab("Unique Genus")
cbPalette2 <- c("deepskyblue", "deepskyblue4")
p <- ggplot(full_dat, aes(x = pneumonia, y = unique_genus)) + geom_boxplot(aes(fill = pneumonia), outlier.shape = NA) + theme_classic()
p + stat_compare_means() + scale_fill_manual(values = c("deepskyblue", "deepskyblue4")) + xlab("Pneumonia") +ylab("Unique Genus") +theme(axis.title.x = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), 

##########
# Testing interactions with FEV1 scores 
#########

new_meta <- read.csv("~/all_patient_df_20220802.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
new_meta <- new_meta %>% select(patient_3pad, percent_predicted_FEV1, smoking_status) %>% as.data.frame()
new_meta <- new_meta %>% na.omit() %>% as.data.frame()

full_dat <- full_dat %>% as.data.frame()
full_dat$smoking <- NA
full_dat$histology <- NA
full_dat$shannons <- NA
full_dat$simpsons <- NA
full_dat$fev1_score <- NA

shannons_diversity <- counts %>% t() %>% as.data.frame()
shannons_diversity <- diversity(shannons_diversity, index = "shannon") %>% as.data.frame()
colnames(shannons_diversity) <- "diversity"
  
simpsons_diversity <- counts %>% t() %>% as.data.frame()
simpsons_diversity <- diversity(simpsons_diversity, index = "simpson") %>% as.data.frame()
colnames(simpsons_diversity) <- "diversity"

for (i in 1:nrow(full_dat)) {
  pat <- full_dat[i, ]
  samp <- pat$samples
  pat <- pat$samples
  pat <- strsplit(pat, split = "_")
  pat <- pat[[1]][2] %>% as.character()
  
  hist <- meta_data %>% filter(Sample == samp)
  hist <- hist$Histology
  smoking <- meta_data %>% filter(Sample == samp)
  smoking <- smoking$Smoking_Status
  shannons <- shannons_diversity[samp, ] %>% as.data.frame()
  shannons <- shannons$.
  simpsons <- simpsons_diversity[samp, ] %>% as.data.frame()
  simpsons <- simpsons$.
  
  full_dat[i, "smoking"] <- smoking
  full_dat[i, "histology"] <- hist
  full_dat[i, "shannons"] <- shannons
  full_dat[i, "simpsons"] <- simpsons
  
  
  if (pat %in% new_meta$patient_3pad) {
  
  new_tmp <- new_meta %>% filter(patient_3pad == pat)
  fev1_score <- new_tmp$percent_predicted_FEV1
  full_dat[i, "fev1_score"] <- fev1_score
  
  } else {
    full_dat[i, "fev1_score"] <- NA
  }
}

full_dat <- full_dat %>% as.data.frame()

text_size = 20

comparisons <- list(c("Current Smoker", "Never Smoked"), c("Never Smoked", "Ex-Smoker"), c("Never Smoked", "Recent Ex-Smoker"))
comparisons <- list(c("LUAD", "LUSC"))
p <- ggplot(full_dat, aes(x = histology, y = shannons, fill = histology)) + geom_boxplot()
p + stat_compare_means(comparisons = comparisons)  + xlab("Smoking Status") + ylab("Shannons Diversity") + theme_classic() + theme(axis.title.x = element_text(size = rel(1.5)), axis.title.y = element_text(size = rel(1.5)),
                                                                                                           legend.text=element_text(size=rel(1.5)), 
                                                                                                           axis.text.y = element_text(face="bold",size=14),
                                                                                                           axis.text.x = element_text(face="bold",size=14))
p <- ggplot(full_dat, aes(x = fev1_score, y = unique_genus)) + geom_point()
p + geom_smooth(method = "lm", se = TRUE)  
