#########
# This script allows to plot the physico-chemical signatures
# Boxplot for each parameter 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(ggsignif)
library(ggrepel)
library(patchwork)
library(tidyverse)


## Data et metadata
Data <- read.csv("02_Data/01_Exp_data/Physico_chemical_parameters.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)

## Layers for boxplot
Colors <- c("orangered", "lightskyblue")

## Boxplot
p1 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Total.Carbon, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Total.Carbon), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Total.Carbon, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Total.Carbon), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Total Carbon (in %)")+
  ylim(0,3)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1

p2 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Total.Nitrogen, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Total.Nitrogen), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Total.Nitrogen, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Total.Nitrogen), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Total Nitrogen (in %)")+
  ylim(0.025,0.3)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2

p3 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = C.N, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = C.N), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = C.N, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = C.N), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Ratio C:N")+
  ylim(10,14)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3

p4 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Available.Phosphorus, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Available.Phosphorus), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Available.Phosphorus, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Available.Phosphorus), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Available Phosphorus (in ppm)")+
  ylim(25,60)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p4

p5 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = CEC, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = CEC), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = CEC, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = CEC), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Cationic exchange capacity (in meq/100g)")+
  ylim(9,22.5)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p5

p6 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Exchange.AI, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Exchange.AI), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Exchange.AI, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Exchange.AI), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Exchange aluminium (in meq/100g)")+
  ylim(0.3,1)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p6

p7 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = pH, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = pH), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = pH, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = pH), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("pH")+
  ylim(4,6)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p7

p8 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Clay, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Clay), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Clay, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Clay), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Clay (in %)")+
  ylim(0,30)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p8

p9 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Silt, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Silt), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Silt, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Silt), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Silt (in %)")+
  ylim(0,60)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p9

p10 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Sand, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Sand), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Sand, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Sand), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Sand (in %)")+
  ylim(30,100)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p10

p11 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Enrichment.index, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Enrichment.index), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Enrichment.index, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Enrichment.index), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Enrichment index")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p11

p12 <- ggplot(Data)+
  geom_boxplot(aes(x = Antagonism, y = Structure.index, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Structure.index), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Structure.index, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Structure.index), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Structure index")+
  theme_bw()+
  ylim(25,125)+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p12

p100 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + plot_layout(guides = "collect")
p100

ggsave(plot = p100, dpi = 1000, device = "svg", width = 7.5, height = 6, filename = "04_Results/Supplementary_Figure_S4a.svg")
# Supplementary_Figure_S4a 

##############################################################################################################"""
#########
# This script allows to plot the Enrichment index depending on the Structure index (Supplementary_Figure_3b)
# These two indices are based on the nematode community (Ferris et al. 2001) 
#########


#### Data ##############################
# Contain the enrichment index and the structure index values for each soil
Indices <- read.csv("02_Data/01_Exp_data/Enrichment_and_structure_indices.csv", sep = ",", header = TRUE)

#### Enrichment / Structure graph 
## Vector containing the shapes for each location
Shape <- c(15,16,17)
Colors <- c("orangered", "lightskyblue")

## Plot the graph
p1 <- ggplot(Indices)+
  geom_point(aes(x = SI, y = EI, shape = Location, color = Antagonism), size = 3)+
  scale_color_manual("Antagonism", values = Colors, guide = "none")+
  scale_shape_manual("Location", values = Shape, guide = "none")+
  geom_text_repel(aes(label = Soils, x = SI, y = EI), size = 4)+
  geom_vline(xintercept = 50, linetype = "dashed")+
  geom_hline(yintercept = 50, linetype = "dashed")+
  xlim(0,100)+
  xlab("Structure index")+
  ylim(0,100)+
  ylab("Enrichment index")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
p1

ggsave(plot = p1, dpi = 1000, device = "svg", width = 6, height = 4, filename = "04_Results/Supplementary_Figure_S4b.svg")

## These results are presented in the Supplementary_Figure_S4b
