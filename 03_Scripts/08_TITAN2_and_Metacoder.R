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

## TITAN - All sites
# Will be merged with the metacoder results
All_sites <- All_sites[["sppmax"]]
write.csv(All_sites, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_all_sites.csv")

## TITAN - Battambang
Battambang <- Battambang[["sppmax"]] # Give the result table
Battambang <- Battambang %>% arrange(-Battambang[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Battambang <- Battambang[grep("TRUE", Battambang[,"filter"] == 2),] # Keep only pure and reliable taxa 
Battambang <- Battambang[,-c(4:5,9:13,17)] # Remove extra columns 
Battambang[,2:3] <- round(Battambang[,2:3], 2)
Battambang[,9] <- round(Battambang[,9], 2)

write.csv(Battambang, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_battambang.csv")

## TITAN - Rovieng
Rovieng <- Rovieng[["sppmax"]] # Give the result table
Rovieng <- Rovieng %>% arrange(-Rovieng[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Rovieng <- Rovieng[grep("TRUE", Bacteria_titan_rovieng[,"filter"] == 2),] # Kees only pure and reliable taxa 
Rovieng <- Rovieng[,-c(4:5,9:13,17)] # Remove extra columns 
Rovieng[,2:3] <- round(Rovieng[,2:3], 2)
Rovieng[,9] <- round(Rovieng[,9], 2)

write.csv(Rovieng, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_rovieng.csv")

## TITAN - Stung Chinit
Stung_Chinit <- Stung_Chinit[["sppmax"]] # Give the result table
Stung_Chinit  <- Stung_Chinit %>% arrange(-Stung_Chinit[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Stung_Chinit  <- Stung_Chinit[grep("TRUE", Stung_Chinit[,"filter"] == 2),] # Keep only pure and reliable taxa 
Stung_Chinit  <- Stung_Chinit[,-c(4:5,9:13,17)] # Remove extra columns 
Stung_Chinit[,2:3] <- round(Stung_Chinit[,2:3], 2)
Stung_Chinit[,9] <- round(Stung_Chinit[,9], 2)

write.csv(Stung_Chinit, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_stung_chinit.csv")

## These results are presented in the Supplementary_Table_S6
## The results tables for pure and reliable bacterial increasers at each site were merged using Microsoft Excel
## A column was added and named "Antagonism against PPNs" where the bacterial taxa with known antagonistic activity to PPNs were identified 

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

## TITAN - Battambang
Battambang <- Battambang[["sppmax"]] # Give the result table
Battambang <- Battambang %>% arrange(-Battambang[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Battambang <- Battambang[grep("TRUE", Battambang[,"filter"] == 2),] # Keep only pure and reliable taxa 
Battambang <- Battambang[,-c(4:5,9:13,17)] # Remove extra columns 
Battambang[,2:3] <- round(Battambang[,2:3], 2)
Battambang[,9] <- round(Battambang[,9], 2)

write.csv(Battambang, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_battambang.csv")

## TITAN - Rovieng
Rovieng <- Rovieng[["sppmax"]] # Give the result table
Rovieng <- Rovieng %>% arrange(-Rovieng[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Rovieng <- Rovieng[grep("TRUE", Bacteria_titan_rovieng[,"filter"] == 2),] # Kees only pure and reliable taxa 
Rovieng <- Rovieng[,-c(4:5,9:13,17)] # Remove extra columns 
Rovieng[,2:3] <- round(Rovieng[,2:3], 2)
Rovieng[,9] <- round(Rovieng[,9], 2)

write.csv(Rovieng, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_rovieng.csv")

## TITAN - Stung Chinit
Stung_Chinit <- Stung_Chinit[["sppmax"]] # Give the result table
Stung_Chinit  <- Stung_Chinit %>% arrange(-Stung_Chinit[,"ienv.cp"], by.group = TRUE) # Arrange results by change point value 
Stung_Chinit  <- Stung_Chinit[grep("TRUE", Stung_Chinit[,"filter"] == 2),] # Keep only pure and reliable taxa 
Stung_Chinit  <- Stung_Chinit[,-c(4:5,9:13,17)] # Remove extra columns 
Stung_Chinit[,2:3] <- round(Stung_Chinit[,2:3], 2)
Stung_Chinit[,9] <- round(Stung_Chinit[,9], 2)

write.csv(Stung_Chinit, file = "Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_stung_chinit.csv")

## These results are presented in the Supplementary_Table_S7
## The results tables for pure and reliable fungal increasers at each site were merged using Microsoft Excel
## A column was added and named "Antagonism against PPNs" where the fungal taxa with known antagonistic activity to PPNs were identified 


#### Metacoder ##############################
#### Data and metadata ####
Bacteria <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Fungi <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table.csv.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Metadata <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Metadata.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)

#### Bacteria ####
## Standardization and percentage transformation
# Add one "fake count" allowing to calculate fold change (HAS/LAS)
# If this is not done, metacoder cannot calculate the fold change because of the division by zero
Bacteria[,8:72] <- Bacteria[,8:72] +1
Bacteria[,8:72] <- tss(Bacteria[,8:72]) # Standardization
Bacteria[,9:73] <- sweep(Bacteria[,8:72], 2, colSums(Bacteria[,8:72]), `/`) * 100 # Genus are now represented as percentage of the sample 

## Changes taxonomy names
# Metacoder need the following taxonomy type = "k__...;p__...;c__... "
Bacteria[,"Kingdom"] <- gsub("^", "k__", Bacteria[,"Kingdom"])
Bacteria[,"Phylum"] <- gsub("^", "p__", Bacteria[,"Phylum"])
Bacteria[,"Class"] <- gsub("^", "c__", Bacteria[,"Class"])
Bacteria[,"Order"] <- gsub("^", "o__", Bacteria[,"Order"])
Bacteria[,"Family"] <- gsub("^", "f__", Bacteria[,"Family"])
Bacteria[,"Genus"] <- gsub("^", "g__", Bacteria[,"Genus"])

Bacteria[,"ID"] <- row.names(Bacteria)
Bacteria <- Bacteria %>% unite(Taxonomy, 2:7, sep = ";", remove = FALSE) # Creates a column "Taxonomy"

# Keep the taxonomy, the genus and the data
Bacteria <- Bacteria[,c(2,8,9:74)]

## Metacoder
# Create a taxmap object
taxmap <- parse_tax_data(Bacteria[,-2], class_cols = "Taxonomy", class_sep = ";")

# Add Metadata in the taxmap environment
# Metacoder will use the first column value as “Treatment 1” to measure enrichment. When sorting metadata by antagonism, the first value is HAS, which allows to calculate the fold change for each taxon by dividing its abundance in HAS by its abundance in LAS
Metadata <- Metadata %>% arrange(Metadata$Antagonism, .by_group = TRUE)
taxmap$data$Metadata <- Metadata

# Calculate the abundances for each taxa
taxmap$data$Abundance <- calc_taxon_abund(taxmap, "tax_data")

# Calculate the mean abundances for each antagonistic group
taxmap$data$Abundance_mean <- calc_group_mean(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Calculate the differences between the two groups 
taxmap$data$Difference_table <- compare_groups(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Performe a Wilcoxon-Mann-Whitney test with a FDR adjustment of the p-value
taxmap <- mutate_obs(taxmap, "Difference_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

# The two data frames will be pooled to know which taxa correspond to the code given by metacoder
Difference_table <- taxmap$data$Difference_table
Taxa <- taxmap$data$tax_data

# Pool of the data
Bacteria <- right_join(Bacteria[,1:2], Taxa[,1:2], by = "Taxonomy")
Results <- right_join(Bacteria[,2:3], Difference_table, by = "taxon_id")
Results <- Results[-grep("TRUE", is.na(Results[,"Genus"])),]

# Keep one row by genus 
Results <- Results %>% distinct(taxon_id, Genus, treatment_1, treatment_2, log2_median_ratio, median_diff, mean_diff, wilcox_p_value)

# Remove the metacoder taxonomy
Results[,"Genus"] <- gsub("^g__", "", Results[,"Genus"])

write.csv(Results, file = "Cambodian_antagonistic_soils/04_Results/Metacoder_Bacteria_HAS_vs_LAS.csv")

#### Fungi ####
## Standardization and percentage transformation
# Add one "fake count" allowing to calculate fold change
# If this is not done, metacoder cannot calculate the fold change because of the division by zero
Fungi[,8:72] <- Fungi[,8:72] +1
Fungi[,8:72] <- tss(Fungi[,8:72]) # Standardization
Fungi[,9:73] <- sweep(Fungi[,8:72], 2, colSums(Fungi[,8:72]), `/`) * 100 # Genus are represented as percentage of the sample 

## Change taxonomy names
# Metacoder need the following taxonomy type = "k__...;p__...;c__... "
Fungi[,"Kingdom"] <- gsub("^", "k__", Fungi[,"Kingdom"])
Fungi[,"Phylum"] <- gsub("^", "p__", Fungi[,"Phylum"])
Fungi[,"Class"] <- gsub("^", "c__", Fungi[,"Class"])
Fungi[,"Order"] <- gsub("^", "o__", Fungi[,"Order"])
Fungi[,"Family"] <- gsub("^", "f__", Fungi[,"Family"])
Fungi[,"Genus"] <- gsub("^", "g__", Fungi[,"Genus"])

Fungi[,"ID"] <- row.names(Fungi)
Fungi <- Fungi %>% unite(Taxonomy, 2:7, sep = ";", remove = FALSE) # Creates a column "Taxonomy"

# Keep the taxonomy, the genus and the data
Fungi <- Fungi[,c(2,8,9:74)]

## Metacoder
# Creates a taxmap object
taxmap <- parse_tax_data(Fungi[,-2], class_cols = "Taxonomy", class_sep = ";")

# Add Metadata in the taxmap environment
# Metacoder will use the first column value as “Treatment 1” to measure enrichment. When sorting metadata by antagonism, the first value is HAS, which allows to calculate the fold change for each taxon by dividing its abundance in HAS by its abundance in LAS
Metadata <- Metadata %>% arrange(Metadata$Antagonism, .by_group = TRUE)
taxmap$data$Metadata <- Metadata

# Calculate the abundances for each taxa
taxmap$data$Abundance <- calc_taxon_abund(taxmap, "tax_data")

# Calculate the mean abundances for each antagonistic group
taxmap$data$Abundance_mean <- calc_group_mean(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Calculate the differences between the two groups 
taxmap$data$Difference_table <- compare_groups(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonis)

# Performe a Wilcoxon-Mann-Whitney test with a FDR adjustment of the p-value
taxmap <- mutate_obs(taxmap, "Difference_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

# The two data frames will be pooled to know which taxa correspond to the code given by metacoder
Difference_table <- taxmap$data$Difference_table
Taxa <- taxmap$data$tax_data

# Pool of the data
Fungi <- right_join(Fungi[,1:2], Taxa[,1:2], by = "Taxonomy")
Results <- right_join(Fungi[,2:3], Difference_table, by = "taxon_id")
Results <- Results[-grep("TRUE", is.na(Results[,"Genus"])),]

# Keep one row by genus
Results <- Results %>% distinct(taxon_id, Genus, treatment_1, treatment_2, log2_median_ratio, median_diff, mean_diff, wilcox_p_value)

# Remove the metacoder taxonomy
Results[,"Genus"] <- gsub("^g__", "", Results[,"Genus"])

write.csv(Results, file = "Cambodian_antagonistic_soils/04_Results/Metacoder_Fungi_HAS_vs_LAS.csv")


#### Satistic tables - All sites ##############################
#### Data ####
# Contains the results from TITAN2
Bacteria_titan <- read.csv("Cambodian_antagonistic_soils/04_Results/TITAN_Bacteria_all_sites.csv", sep = ",", dec = ".", header = TRUE)
Fungi_titan <- read.csv("Cambodian_antagonistic_soils/04_Results/TITAN_Fungi_all_sites.csv.csv", sep = ",", dec = ".", header = TRUE)

# Contains the results from metacoder
Bacteria_metacoder <- read.csv("Cambodian_antagonistic_soils/04_Results/Metacoder_Bacteria_HAS_vs_LAS.csv", sep = ",", dec = ".", header = TRUE)
Fungi_metacoder <- read.csv("Cambodian_antagonistic_soils/04_Results/Metacoder_Bacteria_HAS_vs_LAS", sep = ",", dec = ".", header = TRUE)

#### Bacteria ####
# Arrange TITAN results by the change point value
Bacteria_titan <- Bacteria_titan %>% arrange(-Bacteria_titan[,"ienv.cp"], by.group = TRUE)
# Keep only pure and reliable taxa 
Bacteria_titan <- Bacteria_titan[grep("TRUE", Bacteria_titan[,"filter"] == 2),]
# Remove extra columns
Bacteria_titan <- Bacteria_titan[,-c(4:5,9:13,17)]

Bacteria_titan <- setnames(Bacteria_titan, old = 1, new = "Genus")
Bacteria_metacoder <- Bacteria_metacoder[,-c(1:2,4:5,8)]

# Pool of the data sets
Result_table_bacteria <- left_join(Bacteria_titan, Bacteria_metacoder, by = "Genus")

# Round the results
Result_table_bacteria[,2:3] <- round(Result_table_bacteria[,2:3], 2)
Result_table_bacteria[,9] <- round(Result_table_bacteria[,9], 2)
Result_table_bacteria[,10:11] <- round(Result_table_bacteria[,10:11], 4)
Result_table_bacteria[,12] <- round(Result_table_bacteria[,12], 6)

write.csv(Result_table_bacteria, "Cambodian_antagonistic_soils/04_Results/TITAN_and_Metacoder_Bacteria.csv")

## These results are presented in the Supplementary_Table_S4
## A column was added and named "Antagonism against PPNs", with Microsoft Excel, where the bacterial taxa with known antagonistic activity to PPNs were identified 

#### Fungi ####
# Arrange TITAN2 results by change point value
Fungi_titan <- Fungi_titan %>% arrange(-Fungi_titan[,"ienv.cp"], by.group = TRUE)
# Keep only pure and reliable taxa 
Fungi_titan <- Fungi_titan[grep("TRUE", Fungi_titan[,"filter"] == 2),]
# Remove extra columns 
Fungi_titan <- Fungi_titan[,-c(4:5,9:13,17)]

Fungi_titan <- setnames(Fungi_titan, old = 1, new = "Genus")

Fungi_metacoder <- Fungi_metacoder[,-c(1:2,4:5,8)]

## Pool of the data sets
Result_table_fungi <- left_join(Fungi_titan, Fungi_metacoder, by = "Genus")

# Round the results
Result_table_fungi[,2:3] <- round(Result_table_fungi[,2:3], 2)
Result_table_fungi[,9] <- round(Result_table_fungi[,9], 2)
Result_table_fungi[,10:11] <- round(Result_table_fungi[,10:11], 4)
Result_table_fungi[,12] <- round(Result_table_fungi[,12], 6)

write.csv(Result_table_fungi, "Cambodian_antagonistic_soils/04_Results/TITAN_and_Metacoder_Fungi.csv")

## These results are presented in the Supplementary_Table_S5
## A column was added and named "Antagonism against PPNs", with Microsoft Excel, where the fungal taxa with known antagonistic activity to PPNs were identified 
