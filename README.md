# Predicting regional mutation burden in cancer genomes using chromatin accessibility (CA) and replication timing (RT)

This repository includes source code and processed datasets for the study: 

_Chromatin accessibility of primary human cancers ties regional mutational processes and signatures with tissues of origin_ . 

Oliver Ocsenas and Jüri Reimand (2022) _in revision_.


Tutorials - Jupyter notebooks

* CA2M_RF.ipynb - random forest models of megabase-scale mutation burden, chromatin accessibility and replication timing
* BigWigtoWindowTutorial.ipynb - mapping chromatin signals to megabase-scale windows
* MAFtoWindowTutorial.ipynb - mapping cancer mutations to megabase-scale windows

Tutorials/data - files needed for tutorials

* TumorCA_RT_MBscale.csv.gz - CA and RT tracks for cancer samples, 1-Mbps resolution
* NormalCA_RT_MBscale.csv.gz - CA and RT tracks for normal tissues and cell lines, 1-Mbps resolution
* PCAWG_SNVbinned_MBscale.csv.gz - mutation burden in whole cancer genomes, 1-Mbps resolution
* PCAWG_breastcancer_SNV.MAF.gz - example file of somatic mutations in breast cancer for creating files above
* TCGA_BRCA_ATACSeq_chr1_2.bw - example file of chromatin accessibility in breast cancer for creating files above (chr 1-2 only)

All_code - entire code repository for the project; use on your own responsibility

oliver.ocsenas [@] gmail.com ; juri.reimand [@] utoronto.ca
