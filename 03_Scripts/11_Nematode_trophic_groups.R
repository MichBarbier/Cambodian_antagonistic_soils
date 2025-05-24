########
# This script allows show the difference in nematode composition at their trophic level in each soil
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(data.table)
library(ggsignif)
library(patchwork)
library(tidyverse)


#### Boxplots ##############################
# Contains the number of nematodes from each trophic groups in the samples
Trophic_groups <- read.csv("02_Data/01_Exp_data/Nematodes_trophic_groups.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)

# Create a column with the total number of nematodes in each sample
Trophic_groups[,"Total"] <- apply(Trophic_groups[,1:4], 1, sum)

# Vector with the colors for the boxplots
Colors <- c("orangered", "lightskyblue")

# Boxplot
p1 <- ggplot(Trophic_groups)+
  geom_boxplot(aes(x = Antagonism, y = Total, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Total), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Total, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Total), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Number of nematode larvae")+
  ggtitle("All trophic groups")+
  theme_bw()+
  theme(axis.title.y = element_text(size = 11), axis.text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(Trophic_groups)+
  geom_boxplot(aes(x = Antagonism, y = Bacterivores, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Bacterivores), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Bacterivores, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Bacterivores), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ggtitle("Bacterivores")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 <- ggplot(Trophic_groups)+
  geom_boxplot(aes(x = Antagonism, y = Fungivores, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Fungivores), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Fungivores, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Fungivores), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ggtitle("Fungivores")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p4 <- ggplot(Trophic_groups)+
  geom_boxplot(aes(x = Antagonism, y = Plant_parasites, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Plant_parasites), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Plant_parasites, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Plant_parasites), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ggtitle("Plant-parasites")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p5 <- ggplot(Trophic_groups)+
  geom_boxplot(aes(x = Antagonism, y = Predators_and_omnivores, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Predators_and_omnivores), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Predators_and_omnivores, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Predators_and_omnivores), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ggtitle("Predators and omnivores")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

Model <- "ABCDE"

p10 <- p1 + p2 + p3 + p4 + p5 + plot_layout(design = Model, guides = "collect")

ggsave(plot = p10, dpi = 1000, device = "svg", width = 12, height = 6, filename = "04_Results/Supplementary_Figure_S6.svg")

## These results are presented in the Supplementary_Figure_S6
