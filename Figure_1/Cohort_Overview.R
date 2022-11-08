# Microbial landscape of NSCLC - Cohort Overview Figure 
# Author - Krupa Thakkar 
# Date - 25th September 2022 

# Loading libraries and packages 
library(dplyr)
library(tidyverse)
library(maftools)
library("readxl")
library(ggplot2)
library(stringi)
library(RColorBrewer)
library(gtable)
library(cowplot)
library(readxl)

# Loading data 
absolute_counts <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(absolute_counts) <- absolute_counts[, 1] 
absolute_counts[, 1] <- NULL
absolute_counts <- absolute_counts %>% as.data.frame()
absolute_counts[absolute_counts < 10] <- 0

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

extended_meta <- read_excel("~/Documents/Project 2 -  Microbial landscape of NSCLC/Data/clinical_meta_data_tracerx_tcga.xlsx", sheet = 1)
old_meta <- read.csv("~/Documents/Project 2 -  Microbial landscape of NSCLC/Data/sample_meta_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame() 

taxmat <- readRDS("~/Documents/QIIME2_pipeline/Tx_data/taxmat.rds") %>% as.data.frame()

full_data <- read.delim("~/Documents/Projects/Microbial_TRACERx_P1/Data/RNA_seq/full_genus_matrix.txt", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(full_data) <- full_data[, 1] 
full_data[, 1] <- NULL
full_data <- full_data %>% as.data.frame()
full_data[full_data < 10] <- 0
full_data[is.na(full_data)] <- 0

##### Step 1 : Lollipop plot with total microbial load (ordered by total microbial laod and histology) + creating the large meta_data #####
meta_data <- meta_data%>% filter(!Batch == "2.2")
meta_data <- meta_data %>% filter(Biological == "Tumour")
samples_to_keep <- meta_data$Sample
patients_to_keep <- meta_data$Patient %>% unique()

patients_to_keep <- patients_to_keep %>% as.data.frame()
colnames(patients_to_keep) <- "Name"
rownames(patients_to_keep) <- patients_to_keep$Name
patient_dat <- patients_to_keep
patient_dat$microbial_load <- NA
patient_dat$multi_region <- NA
patient_dat$histology <- NA
patient_dat$smoking_status <- NA
patient_dat$emphysema <- NA 
patient_dat$pneumonia <- NA
patient_dat$chronic_inflammation <- NA
patient_dat$sex<- NA
patient_dat$streptococcus <- NA
patient_dat$staphylococcus <- NA
patient_dat$klebsiella <- NA

absolute_counts <- absolute_counts %>% select(samples_to_keep)

for (i in rownames(patient_dat)) {
  
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Patient == i)
  samples <- tmp_meta$Sample
  
  if(nrow(tmp_meta) > 1) {

    tmp_counts <- absolute_counts %>% select(samples)
    tmp_counts <- tmp_counts %>% t() %>% as.data.frame()
    tmp_counts$total <- rowSums(tmp_counts)
    tmp_counts$name <- "tmp"
    total <- tmp_counts %>% group_by(name) %>% dplyr::summarise(across(everything(), mean)) %>% as.data.frame()
    total <- total$total
    
    patient_dat[i, "microbial_load"] <- total
    patient_dat[i, "multi_region"] <- TRUE
  } else {
    
    tmp_counts <- absolute_counts %>% select(samples)
    colnames(tmp_counts) <- "tmp"
    total <- sum(tmp_counts$tmp)
    
    patient_dat[i, "microbial_load"] <- total
    patient_dat[i, "multi_region"] <- FALSE

  }
}
patient_dat <- patient_dat %>% as.data.frame()

ever_smokers <- c("Recent Ex-Smoker", "Ex-Smoker", "Current Smoker")

for (i in rownames(patient_dat)) {
  
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Patient == i)
  tmp_meta <- tmp_meta[1, ] %>% as.data.frame()
  
  smoking <- tmp_meta$Smoking_Status
  histology <- tmp_meta$Histology
  
  if (smoking %in% ever_smokers){
    patient_dat[i, "smoking_status"] <- "Ever_Smoker"
    patient_dat[i, "histology"] <- histology
  } else {
    patient_dat[i, "smoking_status"] <- "Never_Smoker"
    patient_dat[i, "histology"] <- histology
  }
}
patient_dat <- patient_dat %>% as.data.frame()

for (i in rownames(patient_dat)) {
  
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Patient == i)
  samples <- tmp_meta$Sample
  
  tmp_old <- old_meta[old_meta$Sample %in% samples, ] 
  tmp_old <- tmp_old[1, ] %>% as.data.frame()
  case <- tmp_old$Case
  
  pneumonia <- tmp_old$pneumonia %>% as.character()
  inflam <- tmp_old$Chronic_Inflamm %>% as.character()
  emphysema <- tmp_old$Emphysema %>% as.character()
  
  extended <- extended_meta %>% filter(Hospital_ID == case)
  extended <- extended[1, ] %>% as.data.frame()
  
  sex <- extended$REGSex %>% as.character()
  stage <- extended$REGNSCLCStage %>% as.character()
  
  patient_dat[i, "sex"] <- sex
  patient_dat[i, "stage"] <- stage
  patient_dat[i, "pneumonia"] <- pneumonia
  patient_dat[i, "chronic_inflammation"] <- inflam
  patient_dat[i, "emphysema"] <- emphysema
}
patient_dat <- patient_dat %>% as.data.frame()

##### Step 2 : Adding genus and species level data for the top bacteria relating to pneumonia  #####
for (i in rownames(patient_dat)) {
  
  i <- i %>% as.character()
  
  tmp_meta <- meta_data %>% filter(Patient == i)
  samples <- tmp_meta$Sample
  
  tmp_counts <- absolute_counts %>% select(samples)
  tmp_counts$total <- rowSums(tmp_counts)
  tmp_counts$genus <- rownames(tmp_counts)
  
  klebsiella <- tmp_counts %>% filter(genus == "Klebsiella") %>% as.data.frame()
  klebsiella <- klebsiella$total 
  patient_dat[i, "klebsiella"] <- klebsiella
  
  streptococcus <- tmp_counts %>% filter(genus == "Streptococcus") %>% as.data.frame()
  streptococcus <- streptococcus$total 
  patient_dat[i, "streptococcus"] <- streptococcus
  
  staphylococcus <- tmp_counts %>% filter(genus == "Staphylococcus") %>% as.data.frame()
  staphylococcus <- staphylococcus$total
  patient_dat[i, "staphylococcus"] <- staphylococcus
  
}
patient_dat <- patient_dat %>% as.data.frame()
patient_dat$klebsiella[patient_dat$klebsiella < 10] <- "ABSENT"
patient_dat$streptococcus[patient_dat$streptococcus < 10] <- "ABSENT"
patient_dat$staphylococcus[patient_dat$staphylococcus < 10] <- "ABSENT"

full_data <- full_data %>% t() %>% as.data.frame()
full_data$sample <- NA

for (i in rownames(full_data)) {
  
  i <- i %>% as.character()
  
  tmp <- str_sub(i, end = -7) %>% as.character()
  
  full_data[i, "sample"] <- tmp
  
}
full_data <- full_data %>% as.data.frame()
full_data <- full_data %>% group_by(sample) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
full_data <- full_data %>% filter(sample %in% samples_to_keep)
rownames(full_data) <- full_data[, 1] 
full_data[, 1] <- NULL
full_data <- full_data %>% as.data.frame()
full_data$patient <- NA

for (i in rownames(full_data)) {
  
  i <- i %>% as.character()
  
  tmp <- str_sub(i, end = 8) %>% as.character()
  
  full_data[i, "patient"] <- tmp
  
}
full_data <- full_data %>% as.data.frame()
full_data <- full_data %>% group_by(patient) %>% dplyr::summarise(across(everything(), sum)) %>% as.data.frame()
rownames(full_data) <- full_data[, 1] 
full_data[, 1] <- NULL
full_data <- full_data %>% as.data.frame()

full_data <- full_data %>% t() %>% as.data.frame()
patients_with_rna <- colnames(full_data)
full_data$genus <- NA
full_data$species <- NA
#genus
for (i in rownames(full_data)) {
  full <- i %>% as.character()
  genus_split <- strsplit(i, split = "g__")
  genus_split <- genus_split[[1]][2]
  genus_split_pt2 <- strsplit(genus_split, split = ";s")
  genus <- genus_split_pt2[[1]][1]
  full_data[i, "genus"] <- genus
}
# Species 
for (i in rownames(full_data)) {
  full <- i %>% as.character()
  species_split <- strsplit(i, split = "s__")
  species_split <- species_split[[1]][2]
  species_split_pt2 <- strsplit(species_split, split = ";t")
  species <- species_split_pt2[[1]][1]
  full_data[i, "species"] <- species
}
full_data <- full_data %>% as.data.frame()



#Staphylococcus
for (i in patients_with_rna) {
  
  i <- i %>% as.character()
  
  genus_tmp <- patient_dat[i, ]
  
  genus_tmp <- genus_tmp$staphylococcus
  
  if (genus_tmp == "ABSENT") {
    
  patient_dat[i, "staphylococcus"] <- "ABSENT"
  
  } else {
    tmp_species <- full_data %>% select(i, "genus", "species") %>% as.data.frame()
    tmp_species <- tmp_species %>% filter(genus == "Staphylococcus")
    tmp_species[tmp_species == 0] <- NA
    tmp_species <- tmp_species %>% na.omit()
    
    if (nrow(tmp_species) == 0) {
      
      patient_dat[i, "staphylococcus"] <- "GENUS_LEVEL"
      
    } else {
      
      tmp_species <- tmp_species %>% filter(species == "Staphylococcus_aureus")
      
      if (nrow(tmp_species) == 0) {
        
        patient_dat[i, "staphylococcus"] <- "OTHER_SPECIES"
      } else {
        patient_dat[i, "staphylococcus"] <- "Staphylococcus_aureus"
      }
      
      
    }
    
  }
  
}
patient_dat <- patient_dat %>% as.data.frame()

#Streptococcus
for (i in patients_with_rna) {
  
  i <- i %>% as.character()
  
  genus_tmp <- patient_dat[i, ]
  
  genus_tmp <- genus_tmp$streptococcus
  
  if (genus_tmp == "ABSENT") {
    
    patient_dat[i, "streptococcus"] <- "ABSENT"
    
  } else {
    tmp_species <- full_data %>% select(i, "genus", "species") %>% as.data.frame()
    tmp_species <- tmp_species %>% filter(genus == "Streptococcus")
    tmp_species[tmp_species == 0] <- NA
    tmp_species <- tmp_species %>% na.omit()
    
    if (nrow(tmp_species) == 0) {
      
      patient_dat[i, "streptococcus"] <- "GENUS_LEVEL"
      
    } else {
      
      tmp_species <- tmp_species %>% filter(species == "Streptococcus_pneumoniae")
      
      if (nrow(tmp_species) == 0) {
        
        patient_dat[i, "streptococcus"] <- "OTHER_SPECIES"
      } else {
        patient_dat[i, "streptococcus"] <- "Streptococcus_pneumoniae"
      }
      
      
    }
    
  }
  
}
patient_dat <- patient_dat %>% as.data.frame()

#Klebsiella 
for (i in patients_with_rna) {
  
  i <- i %>% as.character()
  
  genus_tmp <- patient_dat[i, ]
  
  genus_tmp <- genus_tmp$klebsiella
  
  if (genus_tmp == "ABSENT") {
    
    patient_dat[i, "klebsiella"] <- "ABSENT"
    
  } else {
    tmp_species <- full_data %>% select(i, "genus", "species") %>% as.data.frame()
    tmp_species <- tmp_species %>% filter(genus == "Klebsiella")
    tmp_species[tmp_species == 0] <- NA
    tmp_species <- tmp_species %>% na.omit()
    
    if (nrow(tmp_species) == 0) {
      
      patient_dat[i, "klebsiella"] <- "GENUS_LEVEL"
      
    } else {
      
      tmp_species <- tmp_species %>% filter(species == "Klebsiella_pneumoniae")
      
      if (nrow(tmp_species) == 0) {
        
        patient_dat[i, "klebsiella"] <- "OTHER_SPECIES"
      } else {
        patient_dat[i, "klebsiella"] <- "Klebsiella_pneumoniae"
      }
      
      
    }
    
  }
  
}
patient_dat <- patient_dat %>% as.data.frame()

patient_dat["M_LTX916", "klebsiella"] <- "GENUS_LEVEL"
patient_dat["M_LTX916", "streptococcus"] <- "GENUS_LEVEL"
patient_dat["M_LTX916", "staphylococcus"] <- "GENUS_LEVEL"
patient_dat["B_LTX748", "klebsiella"] <- "GENUS_LEVEL"
patient_dat["B_LTX748", "staphylococcus"] <- "GENUS_LEVEL"

patient_dat <- patient_dat %>% as.data.frame()

##### Step 3 : Creating the figure #####
patient_dat$multi_region <- as.factor(patient_dat$multi_region)
histology <- c("NA", "Other", "LUSC", "LUAD")
patient_dat$histology <- factor(patient_dat$histology, levels = histology)
patient_dat$smoking_status <- as.factor(patient_dat$smoking_status)
patient_dat$emphysema <- as.factor(patient_dat$emphysema)
patient_dat$pneumonia <- as.factor(patient_dat$pneumonia)
patient_dat$chronic_inflammation <- as.factor(patient_dat$chronic_inflammation)
patient_dat$sex <- as.factor(patient_dat$sex)
patient_dat$staphylococcus <- as.factor(patient_dat$staphylococcus)
patient_dat$klebsiella <- as.factor(patient_dat$klebsiella)
patient_dat$streptococcus <- as.factor(patient_dat$streptococcus)
patient_dat$stage <- as.factor(patient_dat$stage)
patient_dat_final <- patient_dat

patient_dat <- patient_dat[order(patient_dat$histology, patient_dat$microbial_load, patient_dat$stage, decreasing = T),]  # Kevin wants to swap hist and kraken reads round!
sorted_hist <- sort(unique(patient_dat$histology))
patient_order <- patient_dat$Name %>% as.list()

patient_dat$Name<- factor(patient_dat$Name, levels = patient_order)
#patient_dat$histology <- factor(patient_dat$histology, levels = sorted_hist)
patient_dat$stage <- factor(patient_dat$stage, levels = sort(unique(patient_dat$stage)))

multiregion_cols <- c("aliceblue", "mediumorchid4")
histology_cols <-  c("darkorange", "red4", "darkblue")
stage_cols = c("lightblue", "skyblue","skyblue2","deepskyblue", "dodgerblue")
sex_cols <- c("rosybrown1", "aliceblue")
smoking_cols <- c("dodgerblue", "powderblue")
emphy_cols <- c("white", "lightpink1") 
inflam_cols <- c("white", "lightpink2") 
pneu_cols <- c("white", "lightpink3") 
staphy_cols = c( "white","darkcyan",
                 "violet", "darkslateblue")
kleb_cols = c("white","darkcyan",
               "darkslateblue", "violet")  
strepto_cols = c("white","darkcyan",
                 "violet", "darkslateblue") #FFFF99")
text_size = 20


l_plot <- ggplot(patient_dat, 
                 aes(x=Name, y= microbial_load)) + 
  geom_point() + 
  geom_segment(aes(x=Name, xend=Name, y=0, yend=microbial_load)) +
  ylab("Total microbial load") +
  theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), text = element_text(size = text_size))
l_plot

histology_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = histology)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Histology", values = histology_cols, na.value = "grey") + # "#8C8C8C")+ # '#f0f0f0') +
  ylab("Histology") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
histology_tiles

mr_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = multi_region)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Multiregion", values = multiregion_cols, na.value = "grey") + # "#8C8C8C")+ # '#f0f0f0') +
  ylab("Multiregion") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
mr_tiles

stage_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = stage)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Stage", 
                    values = stage_cols, na.value = "grey") +
  ylab("Stage") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
stage_tiles

sex_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = sex)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Sex",values = sex_cols, na.value = "grey")  +
  ylab("Sex") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
sex_tiles

smoking_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = smoking_status)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Smoking Status", values = smoking_cols, na.value = "grey") + 
  ylab("Smoking Status") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
smoking_tiles

emph_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = emphysema)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Emphysema",values = emphy_cols, na.value = "grey") + 
  ylab("Emphysema") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
emph_tiles

inflam_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = chronic_inflammation)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Inflammation",values = inflam_cols, na.value = "grey") + # "#8C8C8C")+ # '#f0f0f0') +
  ylab("Inflammation") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
inflam_tiles

pneu_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = pneumonia)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Pneumonia",values = pneu_cols, na.value = "grey") + # "#8C8C8C")+ # '#f0f0f0') +
  ylab("Pneumonia") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
pneu_tiles

staphy_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = staphylococcus)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Staphylococcus",values = staphy_cols, na.value = "grey") + # "#8C8C8C")+ # '#f0f0f0') +
  ylab("Staphylococcus") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
staphy_tiles


kleb_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = klebsiella)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Klebsiella",values = kleb_cols, na.value = "grey") + # "#8C8C8C") + # '#f0f0f0') +
  ylab("Klebsiella") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
kleb_tiles

strepto_tiles <- ggplot(patient_dat, aes(x = Name, y = 1, fill = streptococcus)) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = "Streptococcus",values = strepto_cols, na.value = "grey") + # '#f0f0f0') + 
  ylab("Streptococcus") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = text_size))
strepto_tiles


plots <- plot_grid(l_plot,
                   histology_tiles + theme(legend.position = "none"),
                   mr_tiles + theme(legend.position = "none"),
                   stage_tiles + theme(legend.position = "none"),
                   sex_tiles + theme(legend.position = "none"),
                   smoking_tiles + theme(legend.position = "none"),
                   emph_tiles + theme(legend.position = "none"),
                   inflam_tiles + theme(legend.position = "none"), 
                   pneu_tiles + theme(legend.position = "none"),
                   staphy_tiles + theme(legend.position = "none"),
                   kleb_tiles + theme(legend.position = "none"),
                   strepto_tiles + theme(legend.position = "none"), 
                   align = "v", 
                   rel_heights = c(0.5, 0.05, 0.05,0.05,0.05,0.05,0.05, 0.05, 0.05, 0.05, 0.05, 0.05), 
                   axis = "l", ncol = 1)
plots



all_plots_plots <- plot_grid(plots, all_leg,
                             rel_widths = c(8, 0.2),
                             axis = "v",
                             ncol = 2)
all_plots_plots


