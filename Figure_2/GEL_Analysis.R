###############################################################################
####========== Assessing microbiome proportion in lung cancer ============#####
###############################################################################

# Purpose of analysis:
# Determine whether the proportion of non-human reads is a predictor of survival 
# following CPI treatment in lung cancer.

# Secondary endpoint is investigating whether the proportion of non-human reads
# is associated with tumour infiltrating T cell numbers / immune TME
# Author - Ben Simpson

# Libraries
library(survival)
library(survminer)
library(ggplot2)
library(readr)
library(dplyr)

# Read in clinical an genomic variables
# Open mastersheet
Mastersheet_V18 <- read_delim("~/Mastersheet_V18.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) #%>% select(`Participant Id`,os_status,OS_time_treatment,Tumour_type,`Disease Type`,Supertype_B27,Supertype_B07,Coding_TMB,WG_TMB)


# Subset variables of interest
Surv_var <- Mastersheet_V18  %>% select(OS_time_treatment,os_status,Microbiome_proportion,TCRA.fraction.gc.tumour,Disease.Type,Disease.Sub.Type)

# Only lung adeno and squamous carcinomas
Lung_only <- subset(Surv_var, Disease.Type == "LUNG" & (Disease.Sub.Type == "ADENOCARCINOMA" | Disease.Sub.Type == "SQUAMOUS_CELL"))

# Discretise into high and low based on median
Lung_only$Microbial_load <- ifelse(Lung_only$Microbiome_proportion >= median(na.omit(as.numeric(Lung_only$Microbiome_proportion))),"High","Low")

# Make the plot
fit <- survfit(Surv(OS_time_treatment,os_status) ~ Microbial_load, data = Lung_only)
ggsurvplot(
    fit,                     # survfit object with calculated statistics.
      risk.table = T,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      #conf.int = TRUE,         # show confidence intervals for
      # point estimates of survival curves.
      #xlim = c(0,500),         # present narrower X axis, but not affect
      # survival estimates.
      xlab = "Time in months",   # customize X axis label.
      #break.time.by = 100,     # break X axis in time intervals by 500.
      palette = c("deepskyblue2","darkorange"),
      ggtheme = theme_classic(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = TRUE) 

# write to file
pdf(file="Kaplan_Microbial_load_new_color.pdf",width =5 ,height=5)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = T,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # point estimates of survival curves.
  #xlim = c(0,500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  #break.time.by = 100,     # break X axis in time intervals by 500.
  palette = c("deepskyblue2","darkorange"),
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE) 
dev.off()

# Reverse the colours
pdf(file="Kaplan_Microbial_load_new_color_REV.pdf",width =5 ,height=5)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = T,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # point estimates of survival curves.
  #xlim = c(0,500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  #break.time.by = 100,     # break X axis in time intervals by 500.
  palette = c("darkorange","deepskyblue2"),
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE) 
dev.off()

# Scatterplot to compare T cell infiltrate and micro proportion
p3 <- ggplot(Lung_only, aes(as.numeric(Microbiome_proportion), y=TCRA.fraction.gc.tumour)) +
  geom_point(    color="black",
                 fill="deepskyblue2",
                 shape=21,
                 alpha=0.5,
                 size=4,
                 stroke = 0.5) +
  geom_smooth(method=lm , color="#097C9B", fill="#B0EAFA", se=TRUE) + theme_classic() + xlab("Microbial load")+ ylab("T cell infiltrate")
p4 <-p3 + stat_cor(method = "spearman")#, label.x = max(Coding_TMB)/2, label.y = 150,hjust=0.5)
p4 + scale_y_continuous(trans="log2")+ scale_x_continuous(trans="log2")

ggsave("Correlation_TCRA_Microbiome_Lung_SKYBLUE.pdf",
       last_plot(),
       width=15,
       height=15,
       units="cm")

