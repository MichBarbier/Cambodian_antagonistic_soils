#########
# This script allows to :
# - process the sequencing results from the  18S sequencing
# - obtain the rarefaction curves 
# - keep only the fungal taxa according to Naranjo‐Ortizet al. 2022
# - filter the reads (delete reagent contamination, and keep the most seqeunced)
#########


setwd("C:/[[ YOUR PATH]]/Cambodian_antagonistic_soils/01_Processing_seq_results")


#### Packages ##############################

library(data.table)
library(tidyverse)
library(BiocManager)
library(dada2)
library(vegan)


#### DADA2 ##############################

data_files <- "C:/[[ YOUR PATH]]/Cambodian_antagonistic_soils/01_Processing_seq_results/Sequence_18S"
list.files(data_files)

# List the file containing the sequence forwards are named "..._1.fastq"
Seq_fwd <- sort(list.files(data_files, pattern = "_1.fastq", full.names = TRUE))
# List the file containing the sequence reverses are named "..._2.fastq"
Seq_rev <- sort(list.files(data_files, pattern = "_2.fastq", full.names = TRUE))

# Retrieve sample names based on file names 
Sample_names <- sapply(strsplit(basename(Seq_fwd), "_16s"), `[`, 1)

Fwd_filtre <- file.path(data_files, "filtered", paste0(Sample_names, "_F_filt.fastq.gz"))
names(Fwd_filtre) <- Sample_names

Rev_filtre <- file.path(data_files, "filtered", paste0(Sample_names, "_R_filt.fastq.gz"))
names(Rev_filtre) <- Sample_names

# Filter the reads according to their sequencing quality, removes primer length from reads and removes the low-quality base on the 3' side
Seq_trunc <- filterAndTrim(Seq_fwd, Fwd_filtre, Seq_rev, Rev_filtre, truncLen = c(240,240), trimLeft = c(20,18), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)

# Create a model to correct sequencing errors 
Errors_fwd <- learnErrors(Fwd_filtre, multithread = TRUE)
Errors_rev <- learnErrors(Rev_filtre, multithread = TRUE)

# Apply the model to correct the sequencing errors 
Fwd_DADA <- dada(Fwd_filtre, err = Errors_fwd, multithread = TRUE)
Rev_DADA <- dada(Rev_filtre, err = Errors_rev, multithread = TRUE)

# Merge the forward and reverve reads, with a minimal overlap of 12 bases 
Seq_merged <- mergePairs(Fwd_DADA, Fwd_filtre, Rev_DADA, Rev_filtre, minOverlap = 12, verbose = TRUE)

Seq_tab <- makeSequenceTable(Seq_merged)

# Keep the reads with a length ranging between 344 and 430 bases
Seq_tab <- Seq_tab[,nchar(colnames(Seq_tab)) %in% 344:430]

# Remove chimeric sequences 
Seq_tab_nochim <- removeBimeraDenovo(Seq_tab, method = "consensus", multithread = TRUE, verbose = TRUE)

# Assign the taxonomy to the genus level
Taxa <- assignTaxonomy(Seq_tab_nochim, "C:/[[ YOUR PATH]]/Cambodian_antagonistic_soils/01_Processing_seq_results/SILVA/silva_nr99_v138.1_train_set.fa.gz", tryRC = TRUE, multithread = TRUE)

# Create the count table
Seq_tab_analysis <- as.data.frame(t(Seq_tab_nochim))
Seq_tab_analysis$ID <- rownames(Seq_tab_analysis)
Taxa_2 <- as.data.frame(Taxa)
Taxa_2$ID <- rownames(Taxa)
Data_18S <- merge(Taxa_2, Seq_tab_analysis, by = "ID")

write.csv(Data_18S, file = "Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table_unfiltered.csv")


#### Rarefaction curves ##############################

# Create the rarefaction curves, and give the results as a table
Rarecurve_data_18S <- rarecurve(t(Data_18S[,8:77]), tidy = TRUE)

# Create a “Soil” column, which will be used to color replicates of the same soil with the same colors
Rarecurve_data_18S[,"Soil"] <- gsub("_[1-5]$", "", Rarecurve_data_18S[,"Site"])

# Plot the rarefaction curves
p1 <- ggplot(Rarecurve_data_18S)+
  geom_line(aes(x = Sample, y = Species, group = Site, color = Soil))+
  xlab("Number of reads")+
  ylab("Number of ASVs")+
  theme_bw()

ggsave(plot = p1, dpi = 1000, device = "pdf", width = 12, height = 6, filename = "Cambodian_antagonistic_soils/03_Results/Rarefaction_curves_fungi.pdf")


#### Filtering ##############################

# Keep fungal ASVs according to Naranjo‐Ortizet al. 2022
Chytridio <- grep("Chytridiomycota", Data_18S[,"Phylum"])
Basidio <- grep("Basidiomycota", Data_18S[,"Phylum"])
Mucoro <- grep("Mucoromycota", Data_18S[,"Phylum"])
Crypto <- grep("Cryptomycota", Data_18S[,"Phylum"])
Asco <- grep("Ascomycota", Data_18S[,"Phylum"])
Aphel <- grep("Aphelidea", Data_18S[,"Phylum"])
Blastocladio <- grep("Blastocladiomycota", Data_18S[,"Phylum"])
Neocallismatigo <- grep("Neocallimastigomycota", Data_18S[,"Phylum"])

Fungi <- c(Chytridio, Basidio, Mucoro, Crypto, Asco, Aphel, Blastocladio, Neocallismatigo)

Data_18S <- Data_18S[Fungi,]

# Remove the ASVs sequenced in the Control
Ctrl <- grep("Ctrl", colnames(Data_18S))
Conta <- grep("TRUE", apply(Data_18S[,Ctrl], 1, mean) > 0)
Data_16S <- Data_18S[-Conta,-Ctrl]

# Keep the ASVs sequenced at least 5 times are keep
Filtre <- grep("FALSE", apply(Data_18S[,8:72], 1, sum) >= 5)
Data_18S <- Data_18S[-Filtre,]

rows <- c(1:nrow(Data_18S)) # Assign a unique number to the ASVs
rownames(Data_18S) <- rows

# Correct the name of unassigned taxa 
NA_Phyl <- grep("TRUE", is.na(Data_18S[,"Phylum"]))
Data_18S[NA_Phyl,"Phylum"] <- paste0("Unidentified ", Data_18S[NA_Phyl,"Kingdom"])
Data_18S[NA_Phyl,3:7] <- Data_18S[NA_Phyl,"Phylum"]

NA_Class <- grep("TRUE", is.na(Data_18S[,"Class"]))
Data_18S[NA_Class,"Class"] <- paste0("Unidentified ", Data_18S[NA_Class,"Phylum"])
Data_18S[NA_Class,4:7] <- Data_18S[NA_Class,"Class"]

NA_Ord <- grep("TRUE", is.na(Data_18S[,"Order"]))
Data_18S[NA_Ord,"Order"] <- paste0("Unidentified ", Data_18S[NA_Ord,"Class"])
Data_18S[NA_Ord,5:7] <- Data_18S[NA_Ord,"Order"]

NA_Fam <- grep("TRUE", is.na(Data_18S[,"Family"]))
Data_18S[NA_Fam,"Family"] <- paste0("Unidentified ", Data_18S[NA_Fam,"Order"])
Data_18S[NA_Fam,6:7] <- Data_18S[NA_Fam,"Family"]

NA_Gen<-grep("TRUE", is.na(Data_18S[,"Genus"]))
Data_18S[NA_Gen,"Genus"] <- paste0("Unidentified ", Data_18S[NA_Gen,"Family"])
Data_18S[NA_Gen,7] <- Data_18S[NA_Gen,"Genus"]

write.csv(Data_18S, "Cambodian_antagonistic_soils/02_Data/02_Count_tables/Fungal_count_table.csv")
