########
# This script allows to :
# - plot the beta diversity, by using a NMDS and Bray-Curtis distances 
# - show the effect of location, antagonism, and soil physico-chemical parameters on micobiota and nematode communities
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(hilldiv)
library(tidyverse)
library(vegan)


#### Beta diversity ##############################
#### Data ####
# Contains the filtered count tables and the nematode communities
Bacteria <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Fungi <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table.csv.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Nematodes <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Nematode_communities.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)
# Contain the samples information
Metadata <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Metadata.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)

# Add a "Soils" column in the metadata
Metadata[,"Soils"] <- rownames(Metadata)
Metadata[,"Soils"] <- gsub("_.$", "", Metadata[,"Soils"])

# Add a "Labels" column in the metadata which will be use to plot only one soil label on the NMDS
Metadata[,"Labels"] <- rownames(Metadata)
Metadata[,"Labels"] <- gsub("..._._[2-5]$", "", Metadata[,"Labels"])
Metadata[,"Labels"] <- gsub("_1$", "", Metadata[,"Labels"])

#### Bacteria ####
# Data standardization with TSS method
Bacteria[,8:72] <- tss(Bacteria[,8:72])

Bacteria <- as.data.frame(t(Bacteria[,-1:-7])) # Transpose data and remove taxonomy
Bacteria <- Bacteria %>% mutate(.data = Metadata[,c(1,2,12,13)]) # Add samples informations

## NMDS
# Realize the NMDS. Bray-Curtis dissimilarity is used to calculate distances between samples and using 1,000 permutations
NMDS <- metaMDS(Bacteria[,-1:-4], distance = "bray", try = 1000, trymax = 1000)
NMDS_scores <- scores(NMDS)
NMDS_scores_sites <- as.data.frame(NMDS_scores[["sites"]])

# Give the stress value
NMDS$stress
stressplot(NMDS)
Stress_value <- paste0("Stress: ", round(NMDS$stress, 4))

# PERMANOVA
Distances <- vegdist(Bacteria[,-1:-4], method = "bray")

# PERMANOVA on antagonistic categories 
Results_antagonism <- adonis2(Distances~Antagonism, data = Bacteria, permutations = 1000, method = "bray")
# PERMANOVA on location
Results_location <- adonis2(Distances~Location, data = Bacteria, permutations = 1000, method = "bray")

# Get the F-ratio, R2 and P-value from the PERMANOVA tests 
pvalue_permanova_antagonism <- paste0("PERMANOVA antagonism: ", round(Results_antagonism$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_antagonism$`F`[1], 4), "; R2: ",
round(Results_antagonism$`R2`[1], 4))
pvalue_permanova_location <- paste0("PERMANOVA location: ", round(Results_location$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_location$`F`[1], 4), "; R2: ",
round(Results_location$`R2`[1], 4))

# Ordistep
# Scale the physico-chemical parameters
Metadata[,3:11] <- as.data.frame(scale(Metadata[,3:11]))
Metadata <- Metadata[,-c(1,2,12,13)] # Remove columns containing factors 

## RDA
# Allow to perform the ordistep analysis
# Create a model based on physico-chemical parameters explaining statistically the distribution among a RDA
mod0 <- dbrda(Bacteria[,-1:-4] ~ 1, Metadata, distance = "bray") # Model with intercept only
mod1 <- dbrda(Bacteria[,-1:-4] ~ ., Metadata, distance = "bray") # Model with all explanatory variables

## Ordistep
# Ordination of the physico-chemical parameters 
Ordistep_results <- ordistep(mod0, scope = formula(mod1), perm.max = 1000)
Ordistep_results <- Ordistep_results[["anova"]]

## Envfit
# PERMANOVA calculate how the model determined with the ordistep explains the distribution among NMDS
Results_model <- adonis2(Distances~Available.Phosphorus+Total.Nitrogen+Total.Carbon+Exchange.AI+pH+Silt+CEC+Sand+Clay, data = Metadata, permutations = 1000, method = "bray")

# To plot arrows on NMDS
Results_envfit <- envfit(NMDS_scores_sites, Metadata, formula = Available.Phosphorus+Total.Nitrogen+Total.Carbon+Exchange.AI+pH+Silt+CEC+Sand+Clay, method = "bray", permutations = 1000)
Arrows_envfit <- as.data.frame(Results_envfit[["vectors"]][["arrows"]])
Arrows_envfit[,"Labels"] <- row.names(Arrows_envfit)

## Vectors containing colors and shapes of samples for the NMDS
Colors <- c("orangered", "lightskyblue")
Shape <- c(15,16,17)

## Plots the NMDS with arrows
p1 <- ggplot(NMDS_scores_sites, aes(NMDS1,NMDS2))+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Bacteria$Soils), color = "lightgrey", alpha = 0.4, type = "norm", show.legend = FALSE)+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Bacteria$Antagonism, color = Bacteria$Antagonism), alpha = 0.75, size = 0.75, type = "norm", show.legend = FALSE)+
  geom_point(aes(color = Bacteria$Antagonism, shape = Bacteria$Location), size = 3)+
  geom_segment(data = Arrows_envfit, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), color = "darkgrey", alpha = 0.75, size = 0.5, arrow = arrow(angle = 35, length = unit(0.25, "cm")))+
  geom_text_repel(data = Arrows_envfit, aes(label = Labels, x = NMDS1, y = NMDS2), color = "darkgrey")+
  scale_color_manual(Bacteria$Antagonism, values = Colors, guide = "none")+
  scale_shape_manual(Bacteria$Location, values = Shape, guide = "none")+
  geom_text(aes(label = Bacteria$Labels), fontface = "bold", nudge_x = 0.1, nudge_y = 0.1, size = 3.5)+
  ggtitle("Bacterial microbiota")+
  theme_bw()+
  annotate("text", x = 0.25, y = 2.5, label = Stress_value)+
  annotate("text", x = 0.25, y = 2.3, label = pvalue_permanova_antagonism)+
  annotate("text", x = 0.25, y = 2.1, label = pvalue_permanova_location)+
  theme(axis.title.y = element_text(size = 11), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 10))

ggsave(plot = p1, dpi = 1000, device = "svg", width = 12, height = 6, filename = "Cambodian_antagonistic_soils/04_Results/NMDS_Bacteria.svg")
# This figure was modified with Microsoft PowerPoint to remove multiple overlaps 

# To create the statistic table 
write.csv(Results_antagonism, file = "Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_bacteria.csv")
write.csv(Results_location, file = "Cambodian_antagonistic_soils/04_Results/Permanova_location_bacteria.csv")
write.csv(Ordistep_results, file = "Cambodian_antagonistic_soils/04_Results/Ordistep_results_bacteria.csv")
write.csv(Results_model, file = "Cambodian_antagonistic_soils/04_Results/Permanova_model_bacteria.csv")

#### Fungi ####
# Data standardization with TSS method
Fungi[,8:72] <- tss(Fungi[,8:72])

Fungi <- as.data.frame(t(Fungi[,-1:-7])) # Transpose and remove taxonomy
Fungi <- Fungi %>% mutate(.data = Metadata[,c(1,2,12,13)]) # Add metadata

## NMDS
# Realize the NMDS. Bray-Curtis dissimilarity is used to calculate distances between samples and using 1,000 permutations
NMDS <- metaMDS(Fungi[,-1:-4], distance = "bray", try = 1000, trymax = 1000)
NMDS_scores <- scores(NMDS)
NMDS_scores_sites <- as.data.frame(NMDS_scores[["sites"]])

# Give the stress value
NMDS$stress
stressplot(NMDS)
Stress_value <- paste0("Stress: ", round(NMDS$stress, 4))

# PERMANOVA
Distances <- vegdist(Fungi[,-1:-4], method = "bray")

# PERMANOVA on antagonistic categories 
Results_antagonism <- adonis2(Distances~Antagonism, data = Fungi, permutations = 1000, method = "bray")
# PERMANOVA on location
Results_location <- adonis2(Distances~Location, data = Fungi, permutations = 1000, method = "bray")

# Get the F-ratio, R2 and P-value from the PERMANOVA tests 
pvalue_permanova_antagonism <- paste0("PERMANOVA antagonism: ", round(Results_antagonism$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_antagonism$`F`[1], 4), "; R2: ",
round(Results_antagonism$`R2`[1], 4))
pvalue_permanova_location <- paste0("PERMANOVA location: ", round(Results_location$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_location$`F`[1], 4), "; R2: ",
round(Results_location$`R2`[1], 4))

## Ordistep
# Scale the physico-chemical parameters
Metadata[,3:11] <- as.data.frame(scale(Metadata[,3:11]))
Metadata <- Metadata[,-c(1,2,12,13)]

## RDA
# Allow to perform the ordistep analysis
# Create a model based on physico-chemical parameters explaining statistically the distribution among a RDA
mod0 <- dbrda(Fungi[,-1:-4] ~ 1, Metadata, distance = "bray") # Model with intercept only
mod1 <- dbrda(Fungi[,-1:-4] ~ ., Metadata, distance = "bray") # Model with all explanatory variables

## Ordistep
# Ordination of the physico-chemical parameters 
Ordistep_results <- ordistep(mod0, scope = formula(mod1), perm.max = 1000)
Ordistep_results <- Ordistep_results[["anova"]]

## Envfit
# PERMANOVA calculate how the model determined with the ordistep explains the distribution among NMDS
Results_model <- adonis2(Distances~pH+Total.Carbon+Exchange.AI+Total.Nitrogen+Available.Phosphorus+Silt+Clay+CEC+Sand, data = Metadata, permutations = 1000, method = "bray")

# To plot arrows on NMDS
Results_envfit <- envfit(NMDS_scores_sites, Metadata, formula = pH+Total.Carbon+Exchange.AI+Total.Nitrogen+Available.Phosphorus+Silt+Clay+CEC+Sand, method = "bray", permutations = 1000)
Arrows_envfit <- as.data.frame(Results_envfit[["vectors"]][["arrows"]])
Arrows_envfit[,"Labels"] <- row.names(Arrows_envfit)

## Vector containing colors and shapes of samples for the NMDS
Colors <- c("orangered", "lightskyblue")
Shape <- c(15,16,17)

## Plot the NMDS with arrows
# Modified with Microsoft PowerPoint due to overlaps
p2 <- ggplot(NMDS_scores_sites, aes(NMDS1,NMDS2))+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Fungi$Soils), color = "lightgrey", alpha = 0.4, type = "norm", show.legend = FALSE)+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Fungi$Antagonism, color = Fungi$Antagonism), alpha = 0.75, size = 0.75, type = "norm", show.legend = FALSE)+
  geom_point(aes(color = Fungi$Antagonism, shape = Fungi$Location), size = 3)+
  geom_segment(data = Arrows_envfit, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), color = "darkgrey", alpha = 0.75, size = 0.5, arrow = arrow(angle = 35, length = unit(0.25, "cm")))+
  geom_text_repel(data = Arrows_envfit, aes(label = Labels, x = NMDS1, y = NMDS2), color = "darkgrey")+
  scale_color_manual(Fungi$Antagonism, values = Colors, guide = "none")+
  scale_shape_manual(Fungi$Location, values = Shape, guide = "none")+
  geom_text(aes(label = Fungi$Labels), fontface = "bold", nudge_x = 0.1, nudge_y = 0.1, size = 3.5)+
  ggtitle("Fungal microbiota")+
  theme_bw()+
  annotate("text", x = -0.5, y = -2.5, label = Stress_value)+
  annotate("text", x = -0.5, y = -2.75, label = pvalue_permanova_antagonism)+
  annotate("text", x = -0.5, y = -3, label = pvalue_permanova_location)+
  theme(axis.title.y = element_text(size = 11), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 10))

ggsave(plot = p2, dpi = 1000, device = "svg", width = 12, height = 6, filename = "Cambodian_antagonistic_soils/04_Results/NMDS_Fungi.svg")
# This figure was modified with Microsoft PowerPoint to remove multiple overlaps 

# To create the statistic table 
write.csv(Results_antagonism, file = "Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_fungi.csv")
write.csv(Results_location, file = "Cambodian_antagonistic_soils/04_Results/Permanova_location_fungi.csv")
write.csv(Ordistep_results, file = "Cambodian_antagonistic_soils/04_Results/Ordistep_results_fungi.csv")
write.csv(Results_model, file = "Cambodian_antagonistic_soils/04_Results/Permanova_model_fungi.csv")

#### Nematodes ####
Metadata <- Metadata[-c(grep("_[2-5]$", rownames(Metadata))),] # Remove replicates from metabarcoding 
rownames(Metadata) <- gsub("_1$", "", rownames(Metadata))

# Data standardization with TSS method
Nematodes[,4:16] <- tss(Nematodes[,4:16])

Nematodes <- as.data.frame(t(Nematodes[,-1:-3])) # Transpose and remove taxonomy
Nematodes <- Nematodes %>% mutate(.data = Metadata[,c(1,2,12)]) # Add metadata

## NMDS
# Realize the NMDS. Bray-Curtis dissimilarity is used to calculate distances between samples and using 1,000 permutations
NMDS <- metaMDS(Nematodes[,-1:-3], distance = "bray", try = 1000, trymax = 1000)
NMDS_scores <- scores(NMDS)
NMDS_scores_sites <- as.data.frame(NMDS_scores[["sites"]])

# Give the stress value
NMDS$stress
stressplot(NMDS)
Stress_value <- paste0("Stress: ", round(NMDS$stress, 4))

## PERMANOVA
Distances <- vegdist(Nematodes[,-1:-3], method = "bray")

# PERMANOVA on antagonistic categories 
Results_antagonism <- adonis2(Distances~Antagonism, data = Nematodes, permutations = 1000, method = "bray")
# PERMANOVA on antagonistic location
Results_location <- adonis2(Distances~Location, data = Nematodes, permutations = 1000, method = "bray")

# Get the F-ratio, R2 and P-value from the PERMANOVA tests 
pvalue_permanova_antagonism <- paste0("PERMANOVA antagonism: ", round(Results_antagonism$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_antagonism$`F`[1], 4), "; R2: ",
round(Results_antagonism$`R2`[1], 4))
pvalue_permanova_location <- paste0("PERMANOVA location: ", round(Results_location$`Pr(>F)`[1], 4), "; F-ratio: ", round(Results_location$`F`[1], 4), "; R2: ",
round(Results_location$`R2`[1], 4))

## Ordistep
# Scale the physico-chemical parameters
Metadata[,3:11] <- as.data.frame(scale(Metadata[,3:11]))
Metadata <- Metadata[,-c(1,2,12)]

## RDA
# Allow to perform the ordistep analysis
# Create a model based on physico-chemical parameters explaining statistically the distribution among a RDA
mod0 <- dbrda(Nematodes[,-1:-3] ~ 1, Metadata, distance = "bray") # Model with intercept only
mod1 <- dbrda(Nematodes[,-1:-3] ~ ., Metadata, distance = "bray") # Model with all explanatory variables

## Ordistep
# Ordination of the physico-chemical parameters 
Ordistep_results <- ordistep(mod0, scope = formula(mod1), perm.max = 1000)
Ordistep_results <- Ordistep_results[["anova"]]

## Envfit
# PERMANOVA calculate how the model determined with the ordistep explains the distribution among NMDS
Results_model <- adonis2(Distances~Total.Carbon+Clay, data = Metadata, permutations = 1000, method = "bray")

# To plot arrows on NMDS
Results_envfit <- envfit(NMDS_scores_sites, Metadata[,c("Total.Carbon", "Clay")], formula = Total.Carbon+Clay, method = "bray", permutations = 1000)
Arrows_envfit <- as.data.frame(Results_envfit[["vectors"]][["arrows"]])
Arrows_envfit[,"Labels"] <- row.names(Arrows_envfit)

## Vectors containing colors and shapes of samples for the NMDS
Colors <- c("orangered", "lightskyblue")
Shape <- c(15,16,17)

## Plot the NMDS with arrows
p3 <- ggplot(NMDS_scores_sites, aes(NMDS1,NMDS2))+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Nematodes$Antagonism, color = Nematodes$Antagonism), alpha = 0.75, size = 0.75, type = "norm", show.legend = FALSE)+
  geom_point(aes(color = Nematodes$Antagonism, shape = Nematodes$Location), size = 3)+
  geom_segment(data = Arrows_envfit, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), color = "darkgrey", alpha = 0.75, size = 0.5, arrow = arrow(angle = 35, length = unit(0.25, "cm")))+
  geom_text_repel(data = Arrows_envfit, aes(label = Labels, x = NMDS1, y = NMDS2), color = "darkgrey")+
  scale_color_manual(Nematodes$Antagonism, values = Colors, guide = "none")+
  scale_shape_manual(Nematodes$Location, values = Shape, guide = "none")+
  geom_text(aes(label = Nematodes$Soils), fontface = "bold", nudge_x = 0.1, nudge_y = 0.1, size = 3.5)+
  ggtitle("Nematode communities")+
  theme_bw()+
  annotate("text", x = -0.25, y = -0.95, label = Stress_value)+
  annotate("text", x = -0.25, y = -1.1, label = pvalue_permanova_antagonism)+
  annotate("text", x = -0.25, y = -1.25, label = pvalue_permanova_location)+
  theme(axis.title.y = element_text(size = 11), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 10))

ggsave(plot = p3, dpi = 1000, device = "svg", width = 12, height = 6, filename = "Cambodian_antagonistic_soils/04_Results/NMDS_Nematodes.svg")
# Modified with Microsoft PowerPoint due to overlaps

# To create the statistic table 
write.csv(Results_antagonism, file = "Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_nematodes.csv")
write.csv(Results_location, file = "Cambodian_antagonistic_soils/04_Results/Permanova_location_nematodes.csv")
write.csv(Ordistep_results, file = "Cambodian_antagonistic_soils/04_Results/Ordistep_results_nematodes.csv")
write.csv(Results_model, file = "Cambodian_antagonistic_soils/04_Results/Permanova_model_nematodes.csv")


#### Satistic table - PERMANOVA Results ##############################
# Contains PERMANOVA results on the locations
Bacteria_location <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_location_bacteria.csv", sep = ",", dec = ".", header = TRUE)
Fungi_location <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_location_fungi.csv", sep = ",", dec = ".", header = TRUE)
Nematodes_location <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_location_nematodes.csv", sep = ",", dec = ".", header = TRUE)
# Contains PERMANOVA results on the antagonistic activities
Bacteria_antagonism <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_bacteria.csv", sep = ",", dec = ".", header = TRUE)
Fungi_antagonism <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_fungi.csv", sep = ",", dec = ".", header = TRUE)
Nematodes_antagonism <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_antagonism_nematodes.csv", sep = ",", dec = ".", header = TRUE)

# Create the statistic table
Table_statistics <- as.data.frame(matrix(nrow = 6, ncol = 5))
colnames(Table_statistics) <- c("Community", "Factor", "R2", "F-ratio", "p-value")
Table_statistics[,1] <- c("Bacterial microbiota", "Bacterial microbiota", "Fungal microbiota", "Fungal microbiota", "Nematode communities", "Nematode communities")
Table_statistics[,2] <- c("Location", "Antagonism", "Location", "Antagonism", "Location", "Antagonism")

Table_statistics[1,3:5] <- Bacteria_location[1,4:6]
Table_statistics[2,3:5] <- Bacteria_antagonism[1,4:6]
Table_statistics[3,3:5] <- Fungi_location[1,4:6]
Table_statistics[4,3:5] <- Fungi_antagonism[1,4:6]
Table_statistics[5,3:5] <- Nematodes_location[1,4:6]
Table_statistics[6,3:5] <- Nematodes_antagonism[1,4:6]

Table_statistics[,3:4] <- round(Table_statistics[,3:4], 4)
Table_statistics[,5] <- round(Table_statistics[,5], 6)

write.csv(Table_statistics, file = "Cambodian_antagonistic_soils/04_Results/Permanova_results.csv")
# Table was then formated with Microsoft Excel


#### Satistic table - ORDISTEP Results ##############################
# Contains PERMANOVA results on model determined with the Ordistep for the bacterial, fungal and nematodes communities
Bacteria_model <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_model_bacteria.csv", sep = ",", dec = ".", header = TRUE)
Fungi_model <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_model_fungi.csv", sep = ",", dec = ".", header = TRUE)
Nematodes_model <- read.csv("Cambodian_antagonistic_soils/04_Results/Permanova_model_nematodes.csv", sep = ",", dec = ".", header = TRUE)
# Contains Ordistep results for the bacterial, fungal and nematodes communities
Bacteria_ordistep <- read.csv("Cambodian_antagonistic_soils/04_Results/Ordistep_results_bacteria.csv", sep = ",", dec = ".", header = TRUE)
Fungi_ordistep <- read.csv("Cambodian_antagonistic_soils/04_Results/Ordistep_results_fungi.csv", sep = ",", dec = ".", header = TRUE)
Nematodes_ordistep <- read.csv("Cambodian_antagonistic_soils/04_Results/Ordistep_results_nematodes.csv", sep = ",", dec = ".", header = TRUE)

# Create the statistic table
Table_statistics <- as.data.frame(matrix(nrow = 23, ncol = 5))
colnames(Table_statistics) <- c("", "Parameters", "R2", "F-ratio", "p-value")

Table_statistics[1,1:2] <- c("Bacterial microbiota", "Model")
Table_statistics[11,1:2] <- c("Fungal microbiota", "Model")
Table_statistics[21,1:2] <- c("Nematode communities", "Model")

Table_statistics[1,3:5] <- Bacteria_model[1,4:6]
Table_statistics[11,3:5] <- Fungi_model[1,4:6]
Table_statistics[21,3:5] <- Nematodes_model[1,4:6]

Table_statistics[2:10,2] <- Bacteria_ordistep[,1]
Table_statistics[12:20,2] <- Fungi_ordistep[,1]
Table_statistics[22:23,2] <- Nematodes_ordistep[,1]

Table_statistics[2:10,4:5] <- Bacteria_ordistep[,4:5]
Table_statistics[12:20,4:5] <- Fungi_ordistep[,4:5]
Table_statistics[22:23,4:5] <- Nematodes_ordistep[,4:5]

Table_statistics[,3:4] <- round(Table_statistics[,3:4], 4)
Table_statistics[,5] <- round(Table_statistics[,5], 6)

write.csv(Table_statistics, file = "Cambodian_antagonistic_soils/04_Results/Ordistep_and_model_results.csv")
# Table was then formated with Microsoft Excel
