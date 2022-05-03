# Predicting regional mutation burden in cancer genomes using chromatin accessibility (CA) and replication timing (RT)

This repository includes source code, tutorials, and processed datasets for the study: 

_Chromatin accessibility of primary human cancers ties regional mutational processes and signatures with tissues of origin_ . 

Oliver Ocsenas and Jüri Reimand (2022) _in revision_.


Tutorials - Jupyter notebooks

* 1_BigWigtoWindow.ipynb - mapping chromatin signals to megabase-scale windows
* 2_MAFtoWindow.ipynb - mapping cancer mutations to megabase-scale windows
* 3_CA2M_RF.ipynb - random forest models of megabase-scale mutation burden, chromatin accessibility and replication timing
* 4_CA2M_RF_FeatureSelection_Tutorial.ipynb - selecting significant features predicting mutation rates
* 5_CA2M_RF_SHAPscores.ipynb - computing feature importance scores (SHAP)
* 6_CA2M_RF_EnrichedMutations_Tutorial.ipynb - detecting genomic regions with enriched mutations that are not explained by chromatin and replication timing alone

Tutorials/data - files needed for tutorials

* All_CA_RT_100KB_scale.csv.gz - CA and RT tracks for cancer and normal samples, 100-kbps resolution
* All_CA_RT_1MB_scale.csv.gz - CA and RT tracks for cancer and normal samples, 1-Mbps resolution
* NormalCA_RT_MBscale.csv.gz - CA and RT tracks for normal tissues and cell lines, 1-Mbps resolution
* PCAWG_SNVbinned_100KB_scale.csv.gz 
* PCAWG_SNVbinned_MBscale.csv.gz - mutation burden in whole cancer genomes, 1-Mbps resolution
* PCAWG_breastcancer_SNV.MAF.gz - example file of somatic mutations in breast cancer for creating files above
* SHAP_plot.pdf - example plot of feature importance scores (SHAP)
* TCGA_BRCA_ATACSeq_chr1_2.bw - example file of chromatin accessibility in breast cancer for creating files above (chrs 1-2 only)
* TumorCA_RT_MBscale.csv.gz - CA and RT tracks for cancer samples, 1-Mbps resolution

All_code - entire code repository for the project; use on your own responsibility

Contact: oliver.ocsenas [@] gmail.com ; juri.reimand [@] utoronto.ca
