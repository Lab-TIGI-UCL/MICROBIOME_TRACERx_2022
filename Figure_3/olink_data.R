# Microbial landscape of NSCLC - Immunosuppresive environment 
# Author - Krupa Thakkar 
# Date - 24th April 2023 

# Loading libraries 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Loading data
absolute_counts <- read.csv("~/Files/normalised_full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

proteomics_data <- read.csv("~/Files/MISC/BX0344_NPX_annotated_survival_II.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
proteomics_data <- proteomics_data %>% filter(QC_Warning == "PASS") # QC 
proteomics_data <- proteomics_data %>% filter(Timepoint == "Baseline")

##############
# Step 1: Data reformatting
##############
absolute_counts[absolute_counts < 10] <- 0

meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
absolute_counts <- absolute_counts %>% dplyr::select(samples_to_keep)

antiinflammatory <- c("IL10", "IL10RA", "IL10RB")

new_proteomics <- proteomics_data %>% filter(Assay %in% antiinflammatory)

for (i in rownames(new_proteomics)) {
  
  print(i)
  
  olink_id <- i %>% as.character()
  tracerx_id <- proteomics_data %>% filter(TRACERxID == i)
  tracerx_id <- tracerx_id %>% as.data.frame()
  tracerx_id <- tracerx_id[1, ] %>% as.data.frame()
  tracerx_id <- tracerx_id$REGTrialNo %>% as.character()
  
  
  all_samps <- meta_data %>% filter(Patient == tracerx_id)
  all_samps <- all_samps$Sample
  
  counts <- absolute_counts %>% select(all_samps)
  counts <- mean(colSums(counts))
  
  new_proteomics[i, "microbial_load"] <- counts
}

wide_dat <- new_proteomics %>% select(TRACERxID, Assay, NPX) %>% as.data.frame()
tmp_dat <- wide_dat %>%
  group_by(TRACERxID, Assay) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop')
tmp_dat <- tmp_dat %>% as.data.frame()

tmp_dat_wide <- reshape(tmp_dat, idvar = "TRACERxID", timevar = "Assay", direction = "wide")
rownames(tmp_dat_wide) <- tmp_dat_wide$TRACERxID
tmp_dat_wide[,1] <- NULL
all_vars <- colnames(tmp_dat_wide)
tmp_dat_wide_saved <- tmp_dat_wide

for (i in rownames(tmp_dat_wide)) {
  
  print(i)
  
  olink_id <- i %>% as.character()
  tracerx_id <- proteomics_data %>% filter(TRACERxID == i)
  tracerx_id <- tracerx_id %>% as.data.frame()
  tracerx_id <- tracerx_id[1, ] %>% as.data.frame()
  tracerx_id <- tracerx_id$REGTrialNo %>% as.character()
  
  if (tracerx_id %in% meta_data$Patient) {
  
  all_samps <- meta_data %>% filter(Patient == tracerx_id)
  all_samps <- all_samps$Sample
  
  counts <- absolute_counts %>% select(all_samps)
  counts <- mean(colSums(counts))
  
  tmp_dat_wide[i, "microbial_load"] <- counts
  } else {
    tmp_dat_wide[i, "microbial_load"] <- NA
  }
}

tmp_dat_wide <- tmp_dat_wide %>% as.data.frame()

values <- c("NPX.IL10", "NPX.IL10RA", "NPX.IL10RB")

dat <- matrix(nrow = 3, ncol = 1) %>% as.data.frame()
rownames(dat) <- values
colnames(dat) <- "pval"

tmp_dat_wide <- na.omit(tmp_dat_wide)

for (i in rownames(dat)) {
  test <- tmp_dat_wide %>% select(i, microbial_load)
  colnames(test) <- c("a", "b")
  ct <- cor.test(test$a, test$b)
  pval <- ct$p.value
  dat[i, "pval"] <- pval

}

dat$adjsut <- p.adjust(dat$pval)



p <- ggplot(tmp_dat_wide, aes(x = microbial_load, y = NPX.IL10RB)) + geom_point() 
p + geom_smooth(method = "lm", color = "deepskyblue4") + stat_cor(method = "pearson", label.x = 150, label.y = 2.5) + theme_classic()
