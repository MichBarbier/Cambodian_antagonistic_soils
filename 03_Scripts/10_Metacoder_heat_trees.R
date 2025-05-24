########
# This script allows to :
# - measure the enrichment of the increasers in highly antagonistic soils with Metacoder (Foster et al., 2016)
# - plot the heat trees from Figure_4 and Figure_5
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(data.table)
library(ggVennDiagram) # Remove unite() function 
library(hilldiv)
library(metacoder)
library(tidyverse)


#### Bacteria ##############################
#### Heat tree ####
Bacteria <- read.csv("02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Metadata <- read.csv("02_Data/01_Exp_data/Metadata.csv", sep = ";", dec = ".", header = TRUE)
Titan_all_sites <- read.csv("04_Results/TITAN_Bacteria_all_sites.csv", sep = ",", dec = ".", header = TRUE)

Titan_all_sites <- setnames(Titan_all_sites, old = 1, new = "Genus")
Bacteria <- right_join(Bacteria, Titan_all_sites, by = "Genus")
Bacteria <- Bacteria[grep("2", Bacteria[,"filter"]),] # Keep only pure and reliable increasers 
Bacteria <- Bacteria[,-73:-88]

#### Standardization and percentage transformation
Bacteria[,8:72] <- tss(Bacteria[,8:72]) # Standardization
Bacteria[,8:72] <- sweep(Bacteria[,8:72], 2, colSums(Bacteria[,8:72]), `/`) * 100 # Genus are represented as percentage of the sample 

#### Changes taxonomy names
# Metacoder need the following taxonomy type = "k__...;p__...;c__... "
Bacteria[,"Kingdom"] <- gsub("^", "k__", Bacteria[,"Kingdom"])
Bacteria[,"Phylum"] <- gsub("^", "p__", Bacteria[,"Phylum"])
Bacteria[,"Class"] <- gsub("^", "c__", Bacteria[,"Class"])
Bacteria[,"Order"] <- gsub("^", "o__", Bacteria[,"Order"])
Bacteria[,"Family"] <- gsub("^", "f__", Bacteria[,"Family"])
Bacteria[,"Genus"] <- gsub("^", "g__", Bacteria[,"Genus"])

Bacteria[,"ID"] <- row.names(Bacteria)
Bacteria <- Bacteria %>% unite(Taxonomy, 2:7, sep = ";", remove = FALSE) # Creates a column "Taxonomy"

# Keeps the taxonomy, the genus and the data
Bacteria <- Bacteria[,c(2,9:73)]

#### Metacoder
# Creates a taxmap object
taxmap <- parse_tax_data(Bacteria, class_cols = "Taxonomy", class_sep = ";")

# Add Metadata in the taxmap environment
# Metacoder will use the first column value as “Treatment 1” to measure enrichment. When sorting metadata by antagonism, the first value is HAS, which  allows to calculate thefold change for each taxon by dividing the abundance in HAS by the abundance in LAS
Metadata <- Metadata %>% arrange(Metadata$Antagonism, .by_group = TRUE)
taxmap$data$Metadata <- Metadata

# Calculates the abundances for each taxa
taxmap$data$Abundance <- calc_taxon_abund(taxmap, "tax_data")

# Calculates the mean abundances for each antagonistic group
taxmap$data$Abundance_mean <- calc_group_mean(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Calculates the difference between the two groups 
taxmap$data$Difference_table <- compare_groups(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Performes a Mann-Whitney test with a FDR adjustment of the p-value
taxmap <- mutate_obs(taxmap, "Difference_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

p1 <- heat_tree(taxmap,
          node_label = gsub("^[a-z]__", "", taxon_names),
          node_size = n_obs,
          node_label_size = 15,
          edge_label_size = 15,
          node_color = log2_median_ratio,
          node_color_trans = "linear",
          node_color_interval = c(0, 3),
          edge_color_interval = c(0, 3),
          node_color_range = c("lightgrey", "orange", "orangered"),
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Log 2 ratio of median counts
          (HAS / LAS)",
          layout = "re", initial_layout = "re")
p1

ggsave(plot = p1, dpi = 1000, device = "svg", width = 15, height = 15, filename = "Figure_4.svg")

# These results are presented in the Figure_4
## This figure was modified with Microsoft PowerPoint to remove non-statistically enriched increasers from the heat tree
## This figure was modified with Microsoft PowerPoint to colore increasers sequenced in all HAS

#### Venn diagramm ####
# Allow to know in which HAS were sequenced the increasers
Bacteria <- read.csv("02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)

Bacteria <- Bacteria[,-73:-88]
# Pool the ASVs at the genus level
Bacteria <- Bacteria [,-1:-6] # Only keep the genus
Bacteria <- Bacteria %>% pivot_longer(cols = 2:66, names_to = "Soils", values_to = "Count")
Bacteria <- Bacteria %>% group_by(Genus, Soils) %>% summarise(sum(Count)) # Addition of the counts
Bacteria <- Bacteria %>% pivot_wider(id_cols = Genus, names_from = Soils, values_from = 3)
Bacteria <- as.data.frame(Bacteria) # Previous command doesn't give a data.frame
row.names(Bacteria) <- Bacteria[,"Genus"]
Bacteria <- Bacteria[,-1]

Bacteria[,"Bat_1"] <- ifelse(apply(Bacteria[,grep("Bat_1", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Bat_2"] <- ifelse(apply(Bacteria[,grep("Bat_2", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Bat_3"] <- ifelse(apply(Bacteria[,grep("Bat_3", colnames(Bacteria))], 1, sum) > 1, +1, 0)

Bacteria[,"Rov_1"] <- ifelse(apply(Bacteria[,grep("Rov_1", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Rov_2"] <- ifelse(apply(Bacteria[,grep("Rov_2", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Rov_3"] <- ifelse(apply(Bacteria[,grep("Rov_3", colnames(Bacteria))], 1, sum) > 1, +1, 0)

Bacteria[,"Stu_1"] <- ifelse(apply(Bacteria[,grep("Stu_1", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_2"] <- ifelse(apply(Bacteria[,grep("Stu_2", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_3"] <- ifelse(apply(Bacteria[,grep("Stu_3", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_4"] <- ifelse(apply(Bacteria[,grep("Stu_4", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_5"] <- ifelse(apply(Bacteria[,grep("Stu_5", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_6"] <- ifelse(apply(Bacteria[,grep("Stu_6", colnames(Bacteria))], 1, sum) > 1, +1, 0)
Bacteria[,"Stu_7"] <- ifelse(apply(Bacteria[,grep("Stu_7", colnames(Bacteria))], 1, sum) > 1, +1, 0)

Venn <- Bacteria[,66:78]

Rov_1 <- rownames(Venn[grep("TRUE", Venn[,"Rov_1"] == 1),])
Stu_1 <- rownames(Venn[grep("TRUE", Venn[,"Stu_1"] == 1),])
Stu_2 <- rownames(Venn[grep("TRUE", Venn[,"Stu_2"] == 1),])
Stu_3 <- rownames(Venn[grep("TRUE", Venn[,"Stu_3"] == 1),])
Stu_4 <- rownames(Venn[grep("TRUE", Venn[,"Stu_4"] == 1),])
Stu_5 <- rownames(Venn[grep("TRUE", Venn[,"Stu_5"] == 1),])
Stu_6 <- rownames(Venn[grep("TRUE", Venn[,"Stu_6"] == 1),])
Stu_7 <- rownames(Venn[grep("TRUE", Venn[,"Stu_7"] == 1),])

Venn_list <- list(Rov_1, Stu_1, Stu_2, Stu_3, Stu_4, Stu_5, Stu_6, Stu_7)

Venn <- Venn(Venn_list, names = c("Rov_1", "Stu_1", "Stu_2", "Stu_3", "Stu_4", "Stu_5", "Stu_6", "Stu_7"))

Venn_Results <- process_region_data(Venn)
Venn_Results <- Venn_Results %>% arrange(-Venn_Results[,"count"])
Venn_Results <- as.data.frame(Venn_Results[grep("TRUE", Venn_Results["count"] > 0),])
Venn_Results[,"item"] <- as.character(Venn_Results[,"item"])

write.csv(Venn_Results, file = "04_Results/Venn_results_bacteria.csv")



#### Fungi ##############################
#### Heat tree ####
## Data
Fungi <- read.csv("02_Data/02_Count_tables/Fungal_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Metadata <- read.csv("02_Data/01_Exp_data/Metadata.csv.csv", sep = ";", dec = ".", header = TRUE)
Titan_all_sites <- read.csv("04_Results/TITAN_Fungi_all_sites.csv", sep = ",", dec = ".", header = TRUE)

Titan_all_sites <- setnames(Titan_all_sites, old = 1, new = "Genus")
Fungi <- right_join(Fungi, Titan_all_sites, by = "Genus")
Fungi <- Fungi[grep("2", Fungi[,"filter"]),] # Keep only pure and reliable increasers 
Fungi <- Fungi[,-73:-88]

#### Standardization and percentage transformation
Fungi[,8:72] <- tss(Fungi[,8:72]) # Standardization
Fungi[,8:72] <- sweep(Fungi[,8:72], 2, colSums(Fungi[,8:72]), `/`) * 100 # Genus are represented as percentage of the sample 

#### Changes taxonomy names
# Metacoder need the following taxonomy type = "k__...;p__...;c__... "
Fungi[,"Kingdom"] <- gsub("^", "k__", Fungi[,"Kingdom"])
Fungi[,"Phylum"] <- gsub("^", "p__", Fungi[,"Phylum"])
Fungi[,"Class"] <- gsub("^", "c__", Fungi[,"Class"])
Fungi[,"Order"] <- gsub("^", "o__", Fungi[,"Order"])
Fungi[,"Family"] <- gsub("^", "f__", Fungi[,"Family"])
Fungi[,"Genus"] <- gsub("^", "g__", Fungi[,"Genus"])

Fungi[,"ID"] <- row.names(Fungi)
Fungi <- Fungi %>% unite(Taxonomy, 2:7, sep = ";", remove = FALSE) # Creates a column "Taxonomy"

# Keeps the taxonomy, the genus and the data
Fungi <- Fungi[,c(2,9:73)]

#### Metacoder
# Creates a taxmap object
taxmap <- parse_tax_data(Fungi, class_cols = "Taxonomy", class_sep = ";")

# Add Metadata in the taxmap environment
# Metacoder will use the first column value as “Treatment 1” to measure enrichment. When sorting metadata by antagonism, the first value is HAS, which  allows to calculate thefold change for each taxon by dividing the abundance in HAS by the abundance in LAS
Metadata <- Metadata %>% arrange(Metadata$Antagonism, .by_group = TRUE)
taxmap$data$Metadata <- Metadata

# Calculates the abundances for each taxa
taxmap$data$Abundance <- calc_taxon_abund(taxmap, "tax_data")

# Calculates the mean abundances for each antagonistic group
taxmap$data$Abundance_mean <- calc_group_mean(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Calculates the difference between the two groups 
taxmap$data$Difference_table <- compare_groups(taxmap, "Abundance", cols = Metadata$Samples, groups = Metadata$Antagonism)

# Performes a Mann-Whitney test with a FDR adjustment of the p-value
taxmap <- mutate_obs(taxmap, "Difference_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

p1 <- heat_tree(taxmap,
          node_label = gsub("^[a-z]__", "", taxon_names),
          node_size = n_obs,
          node_label_size = 15,
          edge_label_size = 15,
          node_color = log2_median_ratio,
          node_color_trans = "linear",
          node_color_interval = c(0, 2.6),
          edge_color_interval = c(0,2.6),
          node_color_range = c("lightgrey", "orange", "red"),
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Log 2 ratio of median counts
          (HAS / LAS)",
          layout = "re", initial_layout = "re")
p1

ggsave(plot = p1, dpi = 1000, device = "svg", width = 12, height = 6, filename = "Figure_5.svg")

# These results are presented in the Figure_5
## This figure was modified with Microsoft PowerPoint to remove non-statistically enriched increasers from the heat tree
## This figure was modified with Microsoft PowerPoint to colore increasers sequenced in all HAS

#### Venn diagramm ####
# Allow to know in which HAS were sequenced the increasers
Fungi <- read.csv("02_Data/02_Count_tables/Fungal_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)

Fungi <- Fungi[,-73:-88]

# Pool the ASVs at the genus level
Fungi <- Fungi [,-1:-6] # Only keep the genus
Fungi <- Fungi %>% pivot_longer(cols = 2:66, names_to = "Soils", values_to = "Count")
Fungi <- Fungi %>% group_by(Genus, Soils) %>% summarise(sum(Count)) # Addition of the counts
Fungi <- Fungi %>% pivot_wider(id_cols = Genus, names_from = Soils, values_from = 3)
Fungi <- as.data.frame(Fungi) # Previous command doesn't give a data.frame
row.names(Fungi) <- Fungi[,"Genus"]
Fungi <- Fungi[,-1]

Fungi[,"Bat_1"] <- ifelse(apply(Fungi[,grep("Bat_1", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Bat_2"] <- ifelse(apply(Fungi[,grep("Bat_2", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Bat_3"] <- ifelse(apply(Fungi[,grep("Bat_3", colnames(Fungi))], 1, sum) > 1, +1, 0)

Fungi[,"Rov_1"] <- ifelse(apply(Fungi[,grep("Rov_1", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Rov_2"] <- ifelse(apply(Fungi[,grep("Rov_2", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Rov_3"] <- ifelse(apply(Fungi[,grep("Rov_3", colnames(Fungi))], 1, sum) > 1, +1, 0)

Fungi[,"Stu_1"] <- ifelse(apply(Fungi[,grep("Stu_1", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_2"] <- ifelse(apply(Fungi[,grep("Stu_2", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_3"] <- ifelse(apply(Fungi[,grep("Stu_3", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_4"] <- ifelse(apply(Fungi[,grep("Stu_4", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_5"] <- ifelse(apply(Fungi[,grep("Stu_5", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_6"] <- ifelse(apply(Fungi[,grep("Stu_6", colnames(Fungi))], 1, sum) > 1, +1, 0)
Fungi[,"Stu_7"] <- ifelse(apply(Fungi[,grep("Stu_7", colnames(Fungi))], 1, sum) > 1, +1, 0)

Venn <- Fungi[,66:78]

Rov_1 <- rownames(Venn[grep("TRUE", Venn[,"Rov_1"] == 1),])
Stu_1 <- rownames(Venn[grep("TRUE", Venn[,"Stu_1"] == 1),])
Stu_2 <- rownames(Venn[grep("TRUE", Venn[,"Stu_2"] == 1),])
Stu_3 <- rownames(Venn[grep("TRUE", Venn[,"Stu_3"] == 1),])
Stu_4 <- rownames(Venn[grep("TRUE", Venn[,"Stu_4"] == 1),])
Stu_5 <- rownames(Venn[grep("TRUE", Venn[,"Stu_5"] == 1),])
Stu_6 <- rownames(Venn[grep("TRUE", Venn[,"Stu_6"] == 1),])
Stu_7 <- rownames(Venn[grep("TRUE", Venn[,"Stu_7"] == 1),])

Venn_list <- list(Rov_1, Stu_1, Stu_2, Stu_3, Stu_4, Stu_5, Stu_6, Stu_7)

Venn <- Venn(Venn_list, names = c("Rov_1", "Stu_1", "Stu_2", "Stu_3", "Stu_4", "Stu_5", "Stu_6", "Stu_7"))

Venn_Results <- process_region_data(Venn)
Venn_Results <- Venn_Results %>% arrange(-Venn_Results[,"count"])
Venn_Results <- as.data.frame(Venn_Results[grep("TRUE", Venn_Results["count"] > 0),])
Venn_Results[,"item"] <- as.character(Venn_Results[,"item"])

write.csv(Venn_Results, file = "04_Results/Venn_results_fungi.csv")
