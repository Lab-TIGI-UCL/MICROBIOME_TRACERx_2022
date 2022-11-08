# Microbial landscape of NSCLC - Survival volcano plot 
# Author - Krupa Thakkar 
# Date - 29th September 2022 

# Loading packags and libraries 
library(tidyverse)
library(reshape2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(vegan)
library(survminer)
library(survival)
library(cowplot)

# Loading data 
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

survival_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Misc/20220201.TRACERx.all.patients.CTC.release.full.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()

##### Step 1 : Reformatting the data and building a coxph model for microbes assocaiated with reponse ##### 
meta_data <- meta_data %>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_of_interest <- meta_data$Sample
absolute_counts <- absolute_counts %>% select(samples_of_interest)
absolute_counts[absolute_counts < 10] <- 0

new_absolute_counts <- absolute_counts 
new_absolute_counts$genus <- paste0("genus_", 1:nrow(new_absolute_counts))
rownames(new_absolute_counts) <- new_absolute_counts$genus
new_absolute_counts <- new_absolute_counts %>% select(-genus)
new_absolute_counts <- new_absolute_counts %>% t()
new_absolute_counts <- new_absolute_counts %>% as.data.frame()
new_absolute_counts$patient <- NA

for (i in rownames(new_absolute_counts)) {
  i <- i %>% as.character()
  patient <- strsplit(i, split = "_SU")
  patient <- patient[[1]][1] %>% as.character()
  new_absolute_counts[i, "patient"] <- patient
}

new_absolute_counts <- new_absolute_counts %>% group_by(patient) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
rownames(new_absolute_counts) <- new_absolute_counts[, 1] 
new_absolute_counts[, 1] <- NULL
new_absolute_counts <- new_absolute_counts %>% as.data.frame()
new_absolute_counts$patient <- rownames(new_absolute_counts)

survival_pats <- survival_data$REGTrialNo %>% unique()

new_absolute_counts <- new_absolute_counts %>% filter(patient %in% survival_pats)

full_dat <- new_absolute_counts %>% as.data.frame()

for (i in 1:nrow(full_dat)) {
  
  test <- full_dat[i, ] %>% as.data.frame()
  pat <- test$patient %>% as.character()
  
  new_test <- survival_data %>% filter(REGTrialNo == pat)
  ost <- new_test$os_time 
  osc <- new_test$cens_os
  dfst <- new_test$dfs_time
  dfsc <- new_test$cens_dfs
  
  full_dat[i, "OS_time"] <- ost
  full_dat[i, "OS_status"] <- osc
  full_dat[i, "DFS_time"] <- dfst
  full_dat[i, "DFS_status"] <- dfsc
  
}
full_dat <- full_dat %>% as.data.frame()

varlist <- colnames(full_dat[,1:1373])
micro_result <- vector(mode="numeric",length=length(varlist))
HR <- vector(mode="numeric",length=length(varlist))
Upper <- vector(mode="numeric",length=length(varlist))
Lower <- vector(mode="numeric",length=length(varlist))
Name <- vector(mode="numeric",length=length(varlist))

for (i in seq_along(varlist)) {
  
  mod <- as.formula(sprintf("Surv(OS_time, OS_status) ~ %s", varlist[i])) # Sprintf is String print: instead of printing to console, stores output on char buffer which are specified in sprintf. %s represents "insert a string here"
  # generic glm model with imbedded formula
  x <- paste("Model_",i, sep="") # pastes Model, followed by the current iteration value
  assign(x, coxph(formula = mod, data = full_dat)) # assign uses the value x as variable name and assigns the glm output into it
  micro_result[i] <-summary(get(x))$coefficients[1,5]
  HR[i] <-exp(summary(get(x))$coefficients[1,1])
  Upper[i] <-exp(confint(get(x))[1,2])
  Lower[i] <-exp(confint(get(x))[1,1])
  Name[i] <- rownames(summary(get(x))$coefficients)[1]
}

micro_result <- data.frame(Name,HR,Lower,Upper,micro_result)

# Adjust p values and hazard ratio
micro_result$q <- p.adjust(micro_result$micro_result,method="BH")
micro_result$log2HR <- log2(micro_result$HR)
micro_result$genus <- rownames(absolute_counts)
dfs_microresult <- micro_result

# Forest plot 
dfs_microresult <- dfs_microresult %>% na.omit() %>% as.data.frame()
dfs_microresult$prevalence <- NA
rownames(dfs_microresult) <- NULL

new_absolute_counts <- new_absolute_counts %>% as.data.frame()

for (i in 1:nrow(dfs_microresult)) {
  
  new_tmp <- dfs_microresult[i, ] %>% as.data.frame()
  genus <- new_tmp$Name %>% as.character()
  
  new_test <- new_absolute_counts %>% select(genus)
  new_test[new_test < 10] <- NA
  new_test <- new_test %>% na.omit() %>% as.data.frame()
  total <- nrow(new_test)
  prevalence <- (total/210)*100
  dfs_microresult[i, "prevalence"] <- prevalence
  
}

dfs_microresult <- dfs_microresult %>% as.data.frame()
dfs_microresult <- dfs_microresult %>% filter(prevalence > 1)

##### Step 2 : Survival volcano plot ##### 
dfs_microresult <- dfs_microresult %>% as.data.frame()
dfs_microresult$DFS_Associatiion <- NA

p <- ggplot(data=dfs_microresult, aes(x=log2HR, y=-1*log10(q))) + geom_point() + theme_minimal()
p + theme_classic() + xlim(-0.5, 0.5)

EnhancedVolcano(dfs_microresult,
                x = 'log2HR',
                y = 'q', 
                lab = dfs_microresult$genus,
                pCutoff = 0.05, 
                FCcutoff = 0.04, 
                labSize = 3, 
                ylim = c(0, 2), xlim= c(-0.5,0.5), 
                col=c('deeppink', 'deeppink', 'deeppink', 'deepskyblue4'), 
                gridlines.major = FALSE, 
                gridlines.minor = FALSE
                )

