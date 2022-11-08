# Microbial Landscape of NSCLC - Normalisation using Voom-SNM for batch correction 
# Author - Krupa Thakkar 
# Date - 25th September 2022 

# Loading libraries and packages 
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(nlme)
library(plyr)
require(reshape2)
library(gridExtra)
library(cowplot)
library(survival)
library('survminer')
library(scales)
library('dplyr')
library("readr")
library(pheatmap)
library("RColorBrewer")
library(ggplot2)
library(limma)
library(edgeR)
library(gbm)
library(snm)
library(metagenomeSeq)
library(compositions)
library(vegan)

# Loading data 
full_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/full_matrix_decontamination.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
rownames(full_data) <- full_data[,1]
full_data[,1] <- NULL

meta_data <- read.csv("~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Meta_Data/exp_meta.csv", header = TRUE) %>% as.data.frame()

##### Step 1 : Pre-Normalisation PCA #####
samples_meta <- meta_data %>% filter(Type == "Sample") 
samples_meta <- samples_meta %>% filter(!Batch == "2.2")
samples <- samples_meta$Sample
samples_data <- full_data %>% select(samples)
samples_data <- samples_data %>% as.data.frame()

pca_dat <- prcomp(t(samples_data[,2:ncol(samples_data)]))
x<-apply(pca_dat$x, 2, var)
x/sum(x)
pc1_pc2<-data.frame(pca_dat$x[,1:2])
pc1_pc2$Sample<-row.names(pc1_pc2)
pc1_pc2_clin<-merge(pc1_pc2,meta_data,all.x=T,all.y=F)
pc1_pc2_clin$Type <- as.character(pc1_pc2_clin$Type)
pc1_pc2_clin$Batch <- as.character(pc1_pc2_clin$Batch)
# PCA plot
ggplot(pc1_pc2_clin,aes(x=PC1,y=PC2,color= Batch))+geom_point() + theme_classic() + scale_color_manual(values = c("deeppink", "darkorange", "darkgreen"))

##### Step 2 : Post-Normalisation PCA for Voom and SNM #####
counts <- samples_data
dge <- DGEList(counts = counts)
v <- voom(counts, plot=TRUE, normalize="quantile")
v_matrix <- v$E

# PCA post voom
pca_dat<-prcomp(t(v_matrix[,2:ncol(v_matrix)]))
x<-apply(pca_dat$x, 2, var)
x/sum(x)
pc1_pc2<-data.frame(pca_dat$x[,1:2])
pc1_pc2$Sample<-row.names(pc1_pc2)
pc1_pc2_clin<-merge(pc1_pc2,samples_meta,all.x=T,all.y=F)
pc1_pc2_clin$Batch <- as.character(pc1_pc2_clin$Batch)

# PCA plot
ggplot(pc1_pc2_clin,aes(x=PC1,y=PC2,color= Batch))+geom_point()+theme_classic() + stat_ellipse() + scale_color_manual(values = c("deeppink", "darkorange", "darkgreen"))

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

pca_dat<-prcomp(t(snm_data[,2:ncol(snm_data)]))
x<-apply(pca_dat$x, 2, var)
x/sum(x)
pc1_pc2<-data.frame(pca_dat$x[,1:2])
pc1_pc2$Sample<-row.names(pc1_pc2)
#pc1_pc2_clin<-merge(pc1_pc2,meta_data,all.x=T,all.y=F)
pc1_pc2_clin<-merge(pc1_pc2,samples_meta,all.x=T,all.y=F)
pc1_pc2_clin$Batch <- as.character(pc1_pc2_clin$Batch)

# PCA plot
ggplot(pc1_pc2_clin,aes(x=PC1,y=PC2,color= Batch))+geom_point()+theme_classic() +stat_ellipse() + scale_color_manual(values = c("deeppink", "darkorange", "darkgreen"))

##### Step 3 : Saving Voom_SNM matrix #####
#write.csv(snm_data, "~/Documents/Projects/Microbial_TRACERx_P1/Data/16S_seq/Full_Matrix/voom_snm_data.csv")






