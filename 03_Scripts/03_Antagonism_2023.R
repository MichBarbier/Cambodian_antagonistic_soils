#########
# This script allows to plot the boxplot of the antagonism activity in 2023 in each soil
# Antognism : percentage of mobility inhibition 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(multcompView)
library(tidyverse)


#### Antagonism ##############################

# Contain the results of the antagonistic test of 2023 
Antagonism_2023 <- read.csv("02_Data/01_Exp_data/Antagonism_2023.csv", sep = ",", dec = ".", header = TRUE)

# Create a function to determine the statistical groups based on the p-values
tri.to.squ <- function(x)
{
  rn <- row.names(x)
  cn <- colnames(x)
  an <- unique(c(cn,rn))
  myval <-  x[!is.na(x)]
  mymat <-  matrix(1, nrow = length(an), ncol = length(an), dimnames = list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x) == rn[int], colnames(x) == cn[ext]])) next
      mymat[row.names(mymat) == rn[int], colnames(mymat) == cn[ext]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
      mymat[row.names(mymat) == cn[ext], colnames(mymat) == rn[int]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
    }
  }
  return(mymat)
}

# Statistics
# Compare the percentage of mobility between samples with a Wilcoxon-Mann-Whitney test with an adjustment of p-value using a FDR method
Stat_2023 <- pairwise.wilcox.test(Antagonism_2023$Antagonism, Antagonism_2023$Samples, p.adjust.method = "fdr")

# Determine the statistical groups
Mymat <- tri.to.squ(Stat_2023$p.value)

Letters <- multcompLetters(Mymat, compare = "<=", threshold = 0.05, Letters = letters)
Letters <- data.frame(group = names(Letters$Letters), letter = Letters$Letters)

# Vectors containing the colors for the boxplots 
Colors <- c("black", "lightskyblue", "lightskyblue", "lightskyblue", "lightskyblue", "lightskyblue", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered")
Fills <- c("white", "lightskyblue", "lightskyblue", "lightskyblue", "lightskyblue", "lightskyblue", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered", "orangered")

# Plot the boxplots
p1 <- ggplot(Antagonism_2023)+
  geom_boxplot(aes(reorder(Samples, - Antagonism, mean), y = Antagonism, color = reorder(Samples, - Antagonism, mean), fill = reorder(Samples, - Antagonism, mean), alpha = 0.1))+
  scale_fill_manual("Samples", values = Fills, guide = "none")+
  scale_color_manual("Samples", values = Colors, guide = "none")+
  geom_point(aes(x = Samples, y = Antagonism, color = reorder(Samples, - Antagonism, mean)))+
  scale_color_manual("Soil", values = Colors, guide = "none")+
  stat_summary(aes(x = Samples, y = Antagonism), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  geom_text(data = Letters, aes(label = letter, y = 125, x = group), colour = "black", size = 4.5)+
  ylab("Percentage of mobile larvae")+
  scale_y_continuous(breaks = seq(0, 125, 25))+
  xlab("")+
  guides(alpha = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(color = Colors, size = 10), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12))

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 12, height = 6, filename = "04_Results/Figure_1.pdf")

## These results are presented in the Figure_2
