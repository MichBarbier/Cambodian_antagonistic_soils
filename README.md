# Cambodian antagonistic soils
Here is the data and the scripts used in Barbier et al., 2025 : **"Cambodian rice field soils infested by _Meloidogyne_ spp. unveil a potential source of natural pest biocontrol agents"**

**To resume :**
The aims of the study were to identify rice field soils with antagonistic effects against _Meloidogyne graminicola_ and to identfy the microbial taxa responsible for it. Then we sampled soil in 13 rice fields and we characterized their antagonistic effects against the J2 larvae of _M. graminicola_ thanks to an _in vitro_ test where the larvae are challenged with a soil suspension. We also characterized the soil physico-chemical parameters, the nematode communities and sequenced the microbiota focusing on bacteria and fungi. We identified a gradient of antagonistic activity which is lost when the soil suspensions are filtered highlighting a biotic origin. To identify organisms associated with the antagonism we used a Threshold Indicator Taxa ANalysis (Baker and King, 2010) with the TITAN2 package and enrichment analysis with the package metacoder (Foster et al., 2016).

You will find in this repository the data and the scripts allowing to obtain the figures and the tables of the article. But first, you will need to download the raw sequencing data in the ENA database under the accession number PRJEB86942. Once downloaded, register the files ...16s_1.fastq and ...16s_2.fastq in a forlder named **Sequence_16S** and the files ...18s_1.fastq and ...18s_2.fastq in a folder named **Sequences_18S** both placed in the folderd **"01_Processing_seq_results"**. Now you can run the scripts DADA2_Bacteria.R and DADA2_Fungi.R to obtain the count tables which will be used for the different analysis. 

The folder **02_Data** contains the experimental data (metadata, physico-chemical parameters, nematodes communties, antagonism), in 01_Exp_data, and the count tables used for the anlaysis in 02_Count_tables. The folder **03_Scripts** contains all the scripts used to obtain the tables and the figures from the article. The folder **04_Results** is the place where are save the de results (figures and tables) obtained with the scripts. 

The script :
- 01_Sampling_map.R allows to obtain the Supplementary_Figure_S1
- 02_Rarefaction_curves.R allows to obtain the Supplementary_Figure_S2
- 03_Antagonism_2023.R allows to obtain the Figure_1
- 04_Antagonism_2024.R allows to obtain the Supplementary_Figure_S3
- 05_Physico_chemical_parameters.R allows to obtain the Figure_2
- 06_Enrichment_and_structure.R allows to obtain the Supplementary_Figure_S3

####


- 06_Alpha_diversity.R allows to obtain the Supplementary_Figure_S4 and the Supplementary_Table_S2
- 07_Beta_diversity.R allows to obtain the Figure_3, Table_2 and the Supplementary_Table_S3
- 08_TITAN2_and_Metacoder.R allows to obtain the Figure_4, Figure_5, and the Supplementary_Table_S4 to Supplementary_Table_S7
- 09_Nematode_trophic_groups.R allow to obtain the Supplementary_Figure_S5

Some figures have been modified with Microsoft PowerPoint to removes overlaps while some tables have been fomated with Microsoft Excel. 
