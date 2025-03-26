########
# This script allows to :
# - find microbial taxa increasing with the antagonistic gradient with TITAN2 (Baker and King, 2010)
# - measure the enrichment of the increasers in highly antagonistic soils with Metacoder (Foster et al., 2016)
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(data.table)
library(hilldiv)
library(metacoder)
library(tidyverse)
library(TITAN2)


#### TITAN2 ##############################
#### Data ####
# Contains the filtered count tables and the nematode communities
Bacteria <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Fungi <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table.csv.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)

# Contains the antagonistic values as percentage of mobility inhibition 
Antagonism <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Gradient_of_antagonism.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)

#### Bacteria ####
## Standardization with the TSS method
Bacteria[,8:72] <- tss(Bacteria[,8:72])

## Pool the ASVs at the genus level ####
Bacteria <- Bacteria [,-1:-6] # Only keep the genus
Bacteria <- Bacteria %>% pivot_longer(cols = 2:66, names_to = "Soils", values_to = "Count")
Bacteria <- Bacteria %>% group_by(Genus, Soils) %>% summarise(sum(Count)) # Addition of the counts
Bacteria <- Bacteria %>% pivot_wider(id_cols = Genus, names_from = Soils, values_from = 3)
Bacteria <- as.data.frame(Bacteria) # Previous command doesn't give a data.frame
row.names(Bacteria) <- Bacteria[,"Genus"]
Bacteria <- Bacteria[,-1]

##  Split the data depending on the sampling site 
Battambang <- Bacteria[,grep("Bat", colnames(Bacteria))]
Rovieng <- Bacteria[,grep("Rov", colnames(Bacteria))]
Stung_Chinit <- Bacteria[,grep("Stu", colnames(Bacteria))]

Antagonism_battambang <- Antagonism[grep("Bat", rownames(Antagonism)),]
Antagonism_rovieng <- Antagonism[grep("Rov", rownames(Antagonism)),]
Antagonism_stung_chinit <- Antagonism[grep("Stu", rownames(Antagonism)),]

## Filtering 
# Keep the genus sequenced at least in 3 samples
# TITAN2 needs at least 3 occurrences by taxa
Filter <- grep("TRUE",  rowSums(Bacteria[,0:65] > 0) >= 3)
Bacteria <- Bacteria[Filter,]
Bacteria <- as.data.frame(t(Bacteria))

Filter <- grep("TRUE",  rowSums(Battambang[,0:15] > 0) >= 3)
Battambang <- Battambang[Filter,]
Battambang <- as.data.frame(t(Battambang))

Filter <- grep("TRUE",  rowSums(Rovieng[,0:15] > 0) >= 3)
Rovieng <- Rovieng[Filter,]
Rovieng  <- as.data.frame(t(Rovieng))

Filter <- grep("TRUE",  rowSums(Stung_Chinit[,0:35] > 0) >= 3)
Stung_Chinit <- Stung_Chinit[Filter,]
Stung_Chinit  <- as.data.frame(t(Stung_Chinit))

## TITAN 
# Performe the analysis on pooled sampled sites, and splited sampled sites
# 1000 permutations and 500 bootstraps are used
All_sites <- titan(env = Antagonism, txa = Bacteria, numPerm = 1000, nBoot = 500, ncpus = 6)
Battambang <- titan(env = Antagonism_battambang, txa = Battambang, numPerm = 1000, nBoot = 500, ncpus = 6)
Rovieng <- titan(env = Antagonism_rovieng, txa = Rovieng, numPerm = 1000, nBoot = 500, ncpus = 6)
Stung_Chinit <- titan(env = Antagonism_stung_chinit, txa = Stung_Chinit, numPerm = 1000, nBoot = 500, ncpus = 6)

## Plot results on pooled sites
# Shows only the increaser genera with the highest Z-scores
p1 <- plot_taxa_ridges(All_sites, axis.text.y = 7, axis.title.x = 11, z1 = FALSE, xlabel = "Percentage inhibition of larval mobility", xlim = c(0,100), n_ytaxa = 121)

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 8, height = 4, filename = "Cambodian_antagonistic_soils/04_Results/Figure_4.pdf")

## TITAN results 
# Save the results to make the result tables
All_sites <- All_sites[["sppmax"]]
write.csv(All_sites, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_all_sites.csv")
Battambang <- Battambang[["sppmax"]]
write.csv(Battambang, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_battambang.csv")
Rovieng <- Rovieng[["sppmax"]]
write.csv(Rovieng, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_rovieng.csv")
Stung_Chinit <- Stung_Chinit[["sppmax"]]
write.csv(Stung_Chinit, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_stung_chinit.csv")

#### Fungi ####
## Standardization with the TSS method
Fungi[,8:72] <- tss(Fungi[,8:72])

## Pool the ASVs at the genus level
Fungi <- Fungi[,-1:-6] # Only keep the genus
Fungi <- Fungi %>% pivot_longer(cols = 2:66, names_to = "Soils", values_to = "Count")
Fungi <- Fungi %>% group_by(Genus, Soils) %>% summarise(sum(Count)) # Addition of the counts
Fungi <- Fungi %>% pivot_wider(id_cols = Genus, names_from = Soils, values_from = 3)
Fungi <- as.data.frame(Fungi) # Previous command doesn't give a data.frame
row.names(Fungi) <- Fungi[,"Genus"]
Fungi <- Fungi[,-1]

## Split the data depending on the sampling site 
Battambang <- Fungi[,grep("Bat", colnames(Fungi))]
Rovieng <- Fungi[,grep("Rov", colnames(Fungi))]
Stung_Chinit <- Fungi[,grep("Stu", colnames(Fungi))]

Antagonism_battambang <- Antagonism[grep("Bat", rownames(Antagonism)),]
Antagonism_rovieng <- Antagonism[grep("Rov", rownames(Antagonism)),]
Antagonism_stung_chinit <- Antagonism[grep("Stu", rownames(Antagonism)),]

## Filtering 
# Keep the genus sequenced at least in 3 samples
# TITAN2 needs at least 3 occurrences by taxa
Filter <- grep("TRUE", rowSums(Fungi[,0:65] > 0) >= 3)
Fungi <- Fungi[Filter,]
Fungi <- as.data.frame(t(Fungi))

Filter <- grep("TRUE", rowSums(Battambang[,0:15] > 0) >= 3)
Battambang <- Battambang[Filter,]
Battambang <- as.data.frame(t(Battambang))

Filter <- grep("TRUE", rowSums(Rovieng[,0:15] > 0) >= 3)
Rovieng <- Rovieng[Filter,]
Rovieng  <- as.data.frame(t(Rovieng))

Filter <- grep("TRUE", rowSums(Stung_Chinit[,0:35] > 0) >= 3)
Stung_Chinit <- Stung_Chinit[Filter,]
Stung_Chinit  <- as.data.frame(t(Stung_Chinit))

## TITAN 
# Performe the analysis on pooled sampled sites, and splited sampled sites
# 1000 permutations and 500 bootstraps are used
All_sites <- titan(env = Antagonism, txa = Fungi, numPerm = 1000, nBoot = 500, ncpus = 6)
Battambang <- titan(env = Antagonism_battambang, txa = Battambang, numPerm = 1000, nBoot = 500, ncpus = 6)
Rovieng <- titan(env = Antagonism_rovieng, txa = Rovieng, numPerm = 1000, nBoot = 500, ncpus = 6)
Stung_Chinit <- titan(env = Antagonism_stung_chinit, txa = Stung_Chinit, numPerm = 1000, nBoot = 500, ncpus = 6)

## Plot results on pooled sites
# Shows only the increaser genera
p1 <- plot_taxa_ridges(All_sites, axis.text.y = 8, axis.title.x = 11, z1 = FALSE, xlabel = "Percentage inhibition of larval mobility", xlim = c(0,100), n_ytaxa = 244)

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 8, height = 4, filename = "Cambodian_antagonistic_soils/04_Results/Figure_5.pdf")

## TITAN results 
# Save the results to make the result tables
All_sites <- All_sites[["sppmax"]]
write.csv(All_sites, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_all_sites.csv")
Battambang <- Battambang[["sppmax"]]
write.csv(Battambang, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_battambang.csv")
Rovieng <- Rovieng[["sppmax"]]
write.csv(Rovieng, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_rovieng.csv")
Stung_Chinit <- Stung_Chinit[["sppmax"]]
write.csv(Stung_Chinit, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_stung_chinit.csv")

