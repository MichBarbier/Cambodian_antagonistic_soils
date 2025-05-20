#########
# This script allows to plot the soil types of Cambodia according to Crocker 1962
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(tidyverse)
library(sf)


#### Cambodian map ##############################

# Contain .shp files of the Cambodian map and soil types according to Crocker 1962
# .shp files need to be in a folder with .shx files 
# Data obtained from : https://data.opendevelopmentmekong.net/en/dataset/distribution-soil-type-in-cambodia?type=dataset
Soil_types <- st_read("02_Data/01_Exp_data/Shp_files/Soil_types.shp")

# Contain the coordinates of the sampling sites
Sampling_sites <- read.csv("02_Data/01_Exp_data/Coordinates_of_sampling_sites.csv.csv", dec = ".", sep = ",", header = TRUE)

# Vectors containing the soil type names and their colors for ggplot
Soil_names_legend <- c("Mekong and Tonlé sap", "Acid lithosols", "Alluvial lithosols", "Alumisols", "Basic lithosols", "Brown alluvial soils", "Brown hydromorphic soils", "Coastal complex soils", "Cultural hydromorphic soils", "Grey hydromorphic soils", "Lacustrine alluvial soils", "Latosols", "Planosols", "Plinthite podzols", "Plinthitic hydromorphic soils", "Red-yellow podzols", "Regurs soils")
Soil_colors_legend <- c("lightblue1", "antiquewhite1", "antiquewhite2", "antiquewhite3", "antiquewhite4", "rosybrown2", "rosybrown3", "rosybrown4", "tan2", "tan3", "tan4", "indianred2", "indianred3", "indianred4", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4")

# Plot the map with ggplot and geom_sf
p1 <- ggplot()+
  geom_sf(data = Soil_types, aes(fill = Soil_colors_legend))+
  scale_fill_manual("Soil_types", values = Soil_colors_legend, labels = Soil_names_legend, name = "Soil types")+ # Legend
  geom_sf(data = Soil_types[1,], fill = "antiquewhite1")+ # Acid lithosols
  geom_sf(data = Soil_types[2,], fill = "antiquewhite2")+ # Alluvial lithosols
  geom_sf(data = Soil_types[3,], fill = "antiquewhite3")+ # Allumisols
  geom_sf(data = Soil_types[4,], fill = "antiquewhite4")+ # Basic lithosols
  geom_sf(data = Soil_types[5,], fill = "rosybrown2")+ # Brown alluvial soils
  geom_sf(data = Soil_types[6,], fill = "rosybrown3")+ # Brown hydromorphic soils
  geom_sf(data = Soil_types[7,], fill = "rosybrown4")+ # Coastal complex soils
  geom_sf(data = Soil_types[8,], fill = "tan2")+ #Cultural hydromorphic
  geom_sf(data = Soil_types[10,], fill = "tan3")+ # Grey hydromorphic soils
  geom_sf(data = Soil_types[11,], fill = "tan4")+ # Lacustinre alluvial soil
  geom_sf(data = Soil_types[12,], fill = "indianred2")+ # Latosols
  geom_sf(data = Soil_types[13,], fill = "indianred3")+ # Planosols
  geom_sf(data = Soil_types[14,], fill = "indianred4")+ # Plinthite podzols
  geom_sf(data = Soil_types[15,], fill = "lightsteelblue2")+ #Plinthitic hydromorphic soils
  geom_sf(data = Soil_types[16,], fill = "lightsteelblue3")+ # Red-Yellow podzols
  geom_sf(data = Soil_types[17,], fill = "lightsteelblue4")+ # Regur soils
  geom_sf(data = Soil_types[9,], fill = "lightblue1")+ # Mekong and Tonlé Sap
  geom_point(aes(x = Sampling_sites$Longitude, y= Sampling_sites$Latitude), color = "black", size = 1.5, shape = 16)+
  geom_text(aes(x = Sampling_sites$Longitude, y = Sampling_sites$Latitude, label = Sampling_sites$Sites), size = 4, nudge_y = -10000, color = "black")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(legend.justification.right = "left", legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"))

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 8, height = 4, filename = "04_Results/Supplementary_Figure_S1.pdf")

## This map is presented in the Supplementary_Figure_S1
