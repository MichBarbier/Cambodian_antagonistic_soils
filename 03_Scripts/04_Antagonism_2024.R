#########
# This script allows to plot the boxplot of the antagonism activity in 2024 in Bat_3 and Stu_3 soils
# Antognism : percentage of mobility inhibition 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(ggsignif)
library(patchwork)
library(tidyverse)


#### Results 2024
# Contains the result of the antagonism in 2024 for Bat_3 and Stu_3 in percentage of mobility 
Antagonism_2024 <- read.csv("02_Data/Antagonism_2024.csv", sep = ";", dec = ".", header = TRUE)

# Splits data by location
Antagonism_battambang <- Antagonism_2024[grep("Bat", Antagonism_2024_ext[,1]),]
Antagonism_battambang[,1] <- setnames(Antagonism_battambang, old = 1, new = "Suspension")
Antagonism_battambang[,1] <- gsub("Bat_3$", "Unfiltered", Antagonism_battambang[,1])
Antagonism_battambang[,1] <- gsub("Bat_3_Filtered$", "Filtered", Antagonism_battambang[,1])

Antagonism_stung_chinit <- Antagonism_2024[grep("Stu", Antagonism_2024_ext[,1]),]
Antagonism_stung_chinit[,1] <- setnames(Antagonism_stung_chinit, old = 1, new = "Suspension")
Antagonism_stung_chinit[,1] <- gsub("Stu_3$", "Unfiltered", Antagonism_stung_chinit[,1])
Antagonism_stung_chinit[,1] <- gsub("Stu_3_Filtered$", "Filtered", Antagonism_stung_chinit[,1])

#### Boxplots
## Vectors containing the colors for the boxplots
Colors_battambang <- c("black", "lightskyblue")
Fills_battambang <- c("white", "lightskyblue")

Colors_stung_chinit <- c("black", "orangered")
Fills_stung_chinit <- c("white", "orangered" )

## Plots the boxplots
p1 <- ggplot(Antagonism_battambang)+
  geom_boxplot(aes(x = Suspension, y = Antagonism, color = Suspension, fill = Suspension), alpha = 0.65)+
  scale_fill_manual("Suspension", values = Fills_battambang)+
  scale_color_manual("Suspension", values = Colors_battambang)+
  geom_signif(aes(x = Suspension, y = Antagonism), test = wilcox.test, comparisons = list(c("Unfiltered", "Filtered")), map_signif_level = TRUE)+
  geom_point(aes(x = Suspension, y = Antagonism, color = Suspension))+
  scale_color_manual("Soil", values = Colors_battambang, guide = "none")+
  stat_summary(aes(x = Suspension, y = Antagonism), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Percentage of mobile larvae")+
  scale_y_continuous(breaks = seq(0, 155, 25))+
  ylim(0,155)+
  xlab("")+
  ggtitle("Bat_3")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12))
p1

p2 <- ggplot(Antagonism_stung_chinit)+
  geom_boxplot(aes(x = Suspension, y = Antagonism, color = Suspension, fill = Suspension), alpha = 0.65)+
  scale_fill_manual("Suspension", values = Fills_stung_chinit)+
  scale_color_manual("Suspension", values = Colors_stung_chinit)+
  geom_signif(aes(x = Suspension, y = Antagonism), test = wilcox.test, comparisons = list(c("Unfiltered", "Filtered")), map_signif_level = TRUE)+
  geom_point(aes(x = Suspension, y = Antagonism, color = Suspension))+
  scale_color_manual("Soil", values = Colors_stung_chinit, guide = "none")+
  stat_summary(aes(x = Suspension, y = Antagonism), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Percentage of mobile larvae")+
  scale_y_continuous(breaks = seq(0, 155, 25))+
  ylim(0,155)+
  xlab("")+
  ggtitle("Stu_3")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12))
p2

p10 <- p1 + p2

ggsave(plot = p10, dpi = 1000, device = "svg", width = 12, height = 6, filename = "04_Results/Supplementary_Figure_S3.svg")

## These results are presented in the Supplementary_Figure_S3
