#########
# This script allows to :
# - calculate the alpha diversity indices 
# - plot the alpha diversity 
#########


setwd("C:/[[YOUR PATH]]/Cambodian_antagonistic_soils")


#### Packages ##############################

library(ggsignif)
library(microeco)
library(patchwork)
library(tidyverse)


#### Alpha diversity ##############################

# Contain the filtered count tables and the nematode communities
Bacteria <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Bacterial_count_table.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Fungi <- read.csv("Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table.csv.csv", sep = ";", dec = ".", header = TRUE, row.names = 1)
Nematodes <- read.csv("Cambodian_antagonistic_soils/02_Data/01_Exp_data/Nematode_communities.csv", sep = ",", dec = ".", header = TRUE, row.names = 1)

#### Bacteria ####
# Calculate various alpha diversity indices for bacterial microbiota
Bacterial_alpha_diversity <- microtable$new(otu_table = Bacteria[,8:72], tax_table = Bacteria[,2:7], sample_table = as.data.frame(colnames(Bacteria[,8:72])))
Bacterial_alpha_diversity$cal_alphadiv()
Bacterial_alpha_diversity <- Bacterial_alpha_diversity$alpha_diversity

#### Fungi ####
# Calculate various alpha diversity indices for fungal microbiota
Fungal_alpha_diversity <- microtable$new(otu_table = Fungi[,8:72], tax_table = Fungi[,2:7], sample_table = as.data.frame(colnames(Fungi[,8:72])))
Fungal_alpha_diversity$cal_alphadiv()
Fungal_alpha_diversity <- Fungal_alpha_diversity$alpha_diversity

#### Nematodes ####
# Calculate various alpha diversity indices for nematode communities
Nematode_alpha_diversity <- microtable$new(otu_table = Nematodes[,4:16], tax_table = Nematodes[,1:3], sample_table = as.data.frame(colnames(Nematodes[,4:16])))
Nematode_alpha_diversity$cal_alphadiv()
Nematode_alpha_diversity <- Nematode_alpha_diversity$alpha_diversity

# Add soil categories "HAS and LAS" in the result tables 
Bacterial_alpha_diversity[c(grep("Bat", rownames(Bacterial_alpha_diversity)), grep("Rov_[2-3]", rownames(Bacterial_alpha_diversity))),"Antagonism"] <- "Low antagonism"
Fungal_alpha_diversity[c(grep("Bat", rownames(Fungal_alpha_diversity)), grep("Rov_[2-3]", rownames(Fungal_alpha_diversity))),"Antagonism"] <- "Low antagonism"
Nematode_alpha_diversity[c(grep("Bat", rownames(Nematode_alpha_diversity)), grep("Rov_[2-3]", rownames(Nematode_alpha_diversity))),"Antagonism"] <- "Low antagonism"

Bacterial_alpha_diversity[c(grep("Stu", rownames(Bacterial_alpha_diversity)), grep("Rov_1", rownames(Bacterial_alpha_diversity))),"Antagonism"] <- "High antagonism"
Fungal_alpha_diversity[c(grep("Stu", rownames(Fungal_alpha_diversity)), grep("Rov_1", rownames(Fungal_alpha_diversity))),"Antagonism"] <- "High antagonism"
Nematode_alpha_diversity[c(grep("Stu", rownames(Nematode_alpha_diversity)), grep("Rov_1", rownames(Nematode_alpha_diversity))),"Antagonism"] <- "High antagonism"

#### Boxplots #### 
# Vectors containing the colors for boxplots
Colors <- c("orangered", "lightskyblue")

# Plot
p1 <- ggplot(Bacterial_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Observed, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Observed), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Observed, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Observed), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Observed richness")+
  ggtitle("Bacteria")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(Fungal_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Observed, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Observed), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Observed, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Observed), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Observed richness")+
  ggtitle("Fungi")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 <- ggplot(Nematode_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Observed, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Observed), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Observed, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Observed), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Observed richness")+
  ggtitle("Nematodes")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p4 <- ggplot(Bacterial_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Shannon, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Shannon), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Shannon, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Shannon), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Shannon index")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p5 <- ggplot(Fungal_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Shannon, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Shannon), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Shannon, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Shannon), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Shannon index")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p6 <- ggplot(Nematode_alpha_diversity)+
  geom_boxplot(aes(x = Antagonism, y = Shannon, color = Antagonism, fill = Antagonism), alpha = 0.5)+
  geom_signif(aes(x = Antagonism, y = Shannon), test = wilcox.test, comparisons = list(c("High antagonism", "Low antagonism")), map_signif_level = TRUE)+
  scale_color_manual("Antagonism", values = Colors)+
  scale_fill_manual("Antagonism", values = Colors)+
  geom_point(aes(x = Antagonism, y = Shannon, color = Antagonism))+
  scale_color_manual("Antagonism", values = Colors)+
  stat_summary(aes(x = Antagonism, y = Shannon), fun = mean, colour = "black",  geom = "point", shape = 3, size = 2)+
  ylab("Shannon index")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

Model <- "ABC
DEF"

p10 <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(design = Model, guides = "collect")

ggsave(plot = p10, dpi = 1000, device = "pdf", width = 12, height = 6, filename = "Cambodian_antagonistic_soils/04_Results/Supplementary_Figure_S4.pdf")

#### Statistic table #### 
# Creation of the table containing the alpha diversity indices and statistics 
Table_statistics <- as.data.frame(matrix(nrow = 3, ncol = 10))
colnames(Table_statistics) <- c("Observed HAS mean", "Observed HAS sd", "Observed LAS mean", "Observed LAS sd", " Observed p-value", "Shannon HAS mean", "Shannon HAS sd", "Shannon LAS mean", "Shannon LAS sd", "Shannon p-value")
rownames(Table_statistics) <- c("Bacteria", "Fungi", "Nematodes")

# Results
Table_statistics[1,] <- c(mean(Bacterial_alpha_diversity[grep("High", Bacterial_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Bacterial_alpha_diversity[grep("High", Bacterial_alpha_diversity[,"Antagonism"]),"Observed"]), mean(Bacterial_alpha_diversity[grep("Low", Bacterial_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Bacterial_alpha_diversity[grep("Low", Bacterial_alpha_diversity[,"Antagonism"]),"Observed"]), wilcox.test(Bacterial_alpha_diversity[,"Observed"]~Bacterial_alpha_diversity[,"Antagonism"])[["p.value"]], mean(Bacterial_alpha_diversity[grep("High", Bacterial_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Bacterial_alpha_diversity[grep("High", Bacterial_alpha_diversity[,"Antagonism"]),"Shannon"]), mean(Bacterial_alpha_diversity[grep("Low", Bacterial_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Bacterial_alpha_diversity[grep("Low", Bacterial_alpha_diversity[,"Antagonism"]),"Shannon"]), wilcox.test(Bacterial_alpha_diversity[,"Shannon"]~Bacterial_alpha_diversity[,"Antagonism"])[["p.value"]])
Table_statistics[2,] <- c(mean(Fungal_alpha_diversity[grep("High", Fungal_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Fungal_alpha_diversity[grep("High", Fungal_alpha_diversity[,"Antagonism"]),"Observed"]), mean(Fungal_alpha_diversity[grep("Low", Fungal_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Fungal_alpha_diversity[grep("Low", Fungal_alpha_diversity[,"Antagonism"]),"Observed"]), wilcox.test(Fungal_alpha_diversity[,"Observed"]~Fungal_alpha_diversity[,"Antagonism"])[["p.value"]], mean(Fungal_alpha_diversity[grep("High", Fungal_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Fungal_alpha_diversity[grep("High", Fungal_alpha_diversity[,"Antagonism"]),"Shannon"]), mean(Fungal_alpha_diversity[grep("Low", Fungal_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Fungal_alpha_diversity[grep("Low", Fungal_alpha_diversity[,"Antagonism"]),"Shannon"]), wilcox.test(Fungal_alpha_diversity[,"Shannon"]~Fungal_alpha_diversity[,"Antagonism"])[["p.value"]])
Table_statistics[3,] <- c(mean(Nematode_alpha_diversity[grep("High", Nematode_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Nematode_alpha_diversity[grep("High", Nematode_alpha_diversity[,"Antagonism"]),"Observed"]), mean(Nematode_alpha_diversity[grep("Low", Nematode_alpha_diversity[,"Antagonism"]),"Observed"]), sd(Nematode_alpha_diversity[grep("Low", Nematode_alpha_diversity[,"Antagonism"]),"Observed"]), wilcox.test(Nematode_alpha_diversity[,"Observed"]~Nematode_alpha_diversity[,"Antagonism"])[["p.value"]], mean(Nematode_alpha_diversity[grep("High", Nematode_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Nematode_alpha_diversity[grep("High", Nematode_alpha_diversity[,"Antagonism"]),"Shannon"]), mean(Nematode_alpha_diversity[grep("Low", Nematode_alpha_diversity[,"Antagonism"]),"Shannon"]), sd(Nematode_alpha_diversity[grep("Low", Nematode_alpha_diversity[,"Antagonism"]),"Shannon"]), wilcox.test(Nematode_alpha_diversity[,"Shannon"]~Nematode_alpha_diversity[,"Antagonism"])[["p.value"]])

# Round the results
Table_statistics <- round(Table_statistics, 5)
Table_statistics[,1:4] <- round(Table_statistics[,1:4], 2)
Table_statistics[,6:9] <- round(Table_statistics[,6:9], 2)

write.csv(Table_statistics, file = "Cambodian_antagonistic_soils/04_Results/Alhpa_diversity_results.csv")
# Table was formated with with Microsoft Excel
