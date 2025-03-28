#########
# This script allows to plot the Enrichment index depending on the Structure index
# These two indices are based on the nematode community (Ferris et al. 2001) 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(ggrepel)
library(tidyverse)


#### Graph ##############################

# Contain the enrichment index and the structure index values for each soil
Indices <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Enrichment_and_structure_indices.csv", sep = ",", header = TRUE)

#### Enrichment / Structure graph 
## Vector containing the shapes for each location
Shape <- c(15,16,17)

## Plot the graph
p1 <- ggplot(Indices)+
  geom_point(aes(x = SI, y = EI, shape = Location), color = "black", size = 3)+
  scale_shape_manual("Location", values = Shape)+
  geom_text_repel(aes(label = Soils, x = SI, y = EI), size = 4)+
  geom_vline(xintercept = 50, linetype = "dashed")+
  geom_hline(yintercept = 50, linetype = "dashed")+
  xlim(0,100)+
  xlab("Structure index")+
  ylim(0,100)+
  ylab("Enrichment index")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 6, height = 4, filename = "Cambodian_antagonistic_soils/04_Results/Supplementary_Figure_S2.pdf")

## These results are presented in the Supplementary_Figure_S2
