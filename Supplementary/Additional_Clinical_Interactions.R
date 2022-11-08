matrix <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(matrix) <- matrix[, 1] 
matrix[, 1] <- NULL
matrix <- matrix %>% as.data.frame()
counts <- matrix

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/final_meta.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(meta_data) <- meta_data[, 1]
meta_data[, 1] <- NULL
meta_data <- meta_data %>% as.data.frame()

old_meta <- read.csv("~/Documents/Project 2 -  Microbial landscape of NSCLC/Data/sample_meta_data.csv", header = TRUE, check.names = FALSE) %>% as.data.frame() 

meta_data <- meta_data%>% filter(!Batch %in% "2.2")
samples_to_keep <- c("Tumour")
meta_data <- meta_data %>% filter(Biological %in% samples_to_keep)
samples_to_keep <- meta_data$Sample

counts <- counts %>% dplyr::select(samples_to_keep)
counts[counts < 10] <- 0

meta_data$unique_genus <- NA
meta_data$microbial_load <- NA

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp <- counts %>% select(i) %>% as.data.frame()
  
  unique_genus <- tmp %>% as.data.frame()
  unique_genus[unique_genus == 0] <- NA
  unique_genus <- unique_genus %>% na.omit()
  unique_genus <- nrow(unique_genus) 
  meta_data[i, "unique_genus"] <- unique_genus
  
  microbial_load <- tmp %>% as.data.frame()
  microbial_load[microbial_load == 0] <- NA
  microbial_load <- microbial_load %>% na.omit()
  microbial_load <- colSums(microbial_load) 
  meta_data[i, "microbial_load"] <- microbial_load
  
}
meta_data <- meta_data %>% as.data.frame()

meta_data$smoking_full <- NA
ever_smokers <- c("Recent Ex-Smoker", "Ex-Smoker", "Current Smoker")

for (i in samples_to_keep) {
  i <- i %>% as.character()
  tmp <- meta_data[i, ] %>% as.data.frame()
  smoking <- tmp$Smoking_Status %>% as.character()
  
  if (smoking %in% ever_smokers) {
    meta_data[i, "smoking_full"] <- "Ever_Smoker"
  } else {
    meta_data[i, "smoking_full"] <- "Never_Smoker"
  }
}

meta_data <- meta_data %>% as.data.frame()
meta_data <- meta_data %>% filter(!Sample == "B_LTX748_SU_T1-R2")

p <- ggplot(meta_data, aes(x = Histology, y = microbial_load, fill = Histology)) + geom_boxplot()
p + stat_compare_means() + scale_fill_manual(values = c("lightpink","orchid", "grey")) + xlab("Histology") + ylab("Microbial Load") + theme_classic() + theme(axis.title.x = element_text(size = rel(1.5)), 
                                                                                                                                                                axis.title.y = element_text(size = rel(1.5)), 
                                                                                                                                                                legend.text=element_text(size=rel(1.5)), 
                                                                                                                                                                axis.text.y = element_text(face="bold",size=14),
                                                                                                                                                                axis.text.x = element_text(face="bold",size=14))

p
