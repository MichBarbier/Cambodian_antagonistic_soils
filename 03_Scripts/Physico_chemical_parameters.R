#########
# This script allows to plot a Principal Component Analysis of the sampled soils 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(ggrepel)
library(tidyverse)
library(vegan)


#### Principal component analysis ##############################
# Data et metadata
Data_pca <- read.csv("02_Data/01_Exp_data/Physico_chemical_parameters.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)

# Scale the physico-chemical parameters
Data_pca[,3:14] <- scale(Data_pca[,3:14]) # The two first columns are factors

# Makes the PCA
PCA <- pca(Data_pca[,3:14])

PCA_scores <- scores(PCA)
PCA_scores_sites <- as.data.frame(PCA_scores[["sites"]])
PCA_scores_species <- as.data.frame(PCA_scores[["species"]])
PCA_scores_species[,"Var"] <- rownames(PCA_scores_species)

# Get the percentages of variance explained by the first two axes 
PCA_summary <- summary(PCA)
PCA_var <- PCA_summary[["cont"]][["importance"]]
Var1 <- paste0("PC1 - ", round(PCA_var["Proportion Explained", "PC1"]*100, 2), " %")
Var2 <- paste0("PC2 - ", round(PCA_var["Proportion Explained", "PC2"]*100, 2), " %")

# Vector containing the shapes for each location
Shape <- c(15,16,17)

# Plot
p1 <- ggplot(PCA_scores_sites, aes(PC1,PC2))+
  geom_hline(yintercept = 0, color = "black", alpha = 0.7, size = 0.75, linetype = "dotdash")+
  geom_vline(xintercept = 0, color = "black", alpha = 0.7, size = 0.75, linetype = "dotdash")+
  geom_segment(data = PCA_scores_species, aes(x = 0, xend = PC1, y = 0, yend = PC2), color = "grey", arrow = arrow(length = unit(0.25, "cm")), alpha = 0.5)+
  geom_text_repel(data = PCA_scores_species, aes(label = Var, x = PC1, y = PC2), color = "grey", fontface = "bold")+
  geom_point(aes(color = Data_pca$Antagonism, shape = Data_pca$Location), size = 3, color = "black")+
  scale_shape_manual("Location", values = Shape, guide = "none")+
  geom_text_repel(aes(label = rownames(Data_pca)), fontface = "bold", nudge_x = 0.1, nudge_y = 0.09, size = 3)+
  xlab(Var1)+
  ylab(Var2)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

ggsave(plot = p1, dpi = 1000, device = "svg", width = 8, height = 4, filename = "Cambodian_antagonistic_soils/04_Results/Principal_component_analysis.svg")
# This figure was modified with Microsoft PowerPoint to remove multiple overlaps 

## These results are presented in the Figure_1
