#########
# This script allows to plot the rarefaction curves from unfiltrated bacterial and fungal count tables
# This script allows to plot the Supplementary_Figure_2
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(patchwork)
library(tidyverse)
library(vegan)


#### Bacterial rarefaction curves ##############################
# Data
Data_16S <- read.csv("02_Data/02_Count_tables/Bacterial_count_table_unfiltered.csv", sep = ";", header = TRUE, row.names = 1)

#### Rarefaction curve
# Creates the rarefaction curves, and gives the results in a table
Rarecurve_data_16S <- rarecurve(t(Data_16S[,8:77]), tidy = TRUE)

# Create a “Soil” column, which will be used to color replicates of the same soil with the same colors
Rarecurve_data_16S[,"Soil"] <- gsub("_[1-5]$", "", Rarecurve_data_16S[,"Site"])
# Change the name of the control 
Rarecurve_data_16S[,"Soil"] <- gsub("CtrlExtract", " Extraction kit control", Rarecurve_data_16S[,"Soil"])

# Plots the rarefaction curves
p1 <- ggplot(Rarecurve_data_16S)+
  geom_line(aes(x = Sample, y = Species, group = Site, color = Soil))+
  xlab("Number of reads")+
  xlim(0,55000)+
  ylab("Number of ASVs")+
  scale_y_continuous(breaks = seq(0, 2500, 500))+
  ggtitle("Bacterial microbiota")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p1


#### Fungal rarefaction curves ##############################
# Data
Data_18S <- read.csv("02_Data/02_Count_tables/Fungal_count_table_unfiltered.csv", sep = ";", header = TRUE, row.names = 1)

#### Rarefaction curve
# Creates the rarefaction curves, and gives the results in a table
Rarecurve_data_18S <- rarecurve(t(Data_18S[,8:77]), tidy = TRUE)

# Create a “Soil” column, which will be used to color replicates of the same soil with the same colors
Rarecurve_data_18S[,"Soil"] <- gsub("_[1-5]$", "", Rarecurve_data_18S[,"Site"])
# Change the name of the control 
Rarecurve_data_18S[,"Soil"] <- gsub("CtrlExtract", " Extraction kit control", Rarecurve_data_18S[,"Soil"])

# Plots the rarefaction curves
p2 <- ggplot(Rarecurve_data_18S)+
  geom_line(aes(x = Sample, y = Species, group = Site, color = Soil))+
  xlab("Number of reads")+
  xlim(0,55000)+
  ggtitle("Fungal microbiota")+
  ylab("Number of ASVs")+
  theme_bw()
p2


#### Supplementary Figure 2
p10 <- p1 / p2 + plot_layout(guides = "collect")

ggsave(plot = p10, dpi = 1000, device = "svg", width = 12, height = 6, filename = "04_Results/Supplementary_Figure_S2.svg")
