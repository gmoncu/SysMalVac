# SysMalVac
This repository contains the scripts used in the Moncunill _et al._ (2020) manuscript about transcriptional signatures of protection upon malaria sporozoite and RTS,S/AS01E immunizations

The scripts are divided in 4 different blocs of analysis: 

## 1. Differential gene expression analyses

The code is structured in 3 different folders:
- **_CPS_** folder:
  * [`CPS_analysis1.R`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/CPS/CPS_analysis1.R) performs the microarray normalization and filter.
  * [`CPS_analysis2.R`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/CPS/CPS_analysis2.R) performs exploratory differential gene expression not reported in the manuscript. 
  * [`CPS_analysis3.R`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/CPS/CPS_analysis3.R) performes the differential gene expression analysis reported in the manuscript and used for GSEA (point 3 below) and radar plots (point 2 below).
  * [`phenoDATA.txt`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/CPS/phenoDATA.txt) cointains variables from the clinical trial.
- **_RTS,S_** folder: 
  * [`RTSS_analysis1.R`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/RTS%2CS/RTSS_analysis1.R) performs the microarray normalization and filter, and the differential gene expression analysis reported in the manuscript and used for GSEA (point 3 below). 
  * [`RTSS_analysis2.R`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/RTS%2CS/RTSS_analysis2.R) perfomrs the differential gene expression analysis used for radar plots (point 2 below).
  * [`phenoDATA.txt`](https://github.com/gmoncu/SysMalVac/blob/master/1_Differential%20gene%20expression/RTS%2CS/phenoDATA.txt) cointains variables from the clinical trial.
- **_Utils_** folder:
Contains functions used for microarray QC, filter, plots and statistics used in the code of above.

Microarrays data is available at the Gene Expression Omnibus ([GEO](www.ncbi.nlm.nih.gov/geo)) database  (accession no. XXXX).


## 2. PCA, heatmaps, gene correlations and radar charts

   * [`PCA_analysis&figures.R`](https://github.com/gmoncu/SysMalVac/blob/master/2_PCA%2C%20heatmaps%2C%20gene%20correlations%20and%20radar%20charts/PCA_analysis%26figures.R) is the code to perform PCA analysis and manuscript figures.
   * [`Heatmap_figures.R`](https://github.com/gmoncu/SysMalVac/blob/master/2_PCA%2C%20heatmaps%2C%20gene%20correlations%20and%20radar%20charts/Heatmaps_figures.R) is the code to perform the heatmap manuscript figures.
   * [`Gene_correlations.R`](https://github.com/gmoncu/SysMalVac/blob/master/2_PCA%2C%20heatmaps%2C%20gene%20correlations%20and%20radar%20charts/Gene_correlations.R) is the code to perform the analysis and manuscript figures of FC PfRBC/uRBC gene expression correlations with parasitemia and prepatency in the CPS study.
   * [`RadarCharts_FCgeneExpression.R`](https://github.com/gmoncu/SysMalVac/blob/master/2_PCA%2C%20heatmaps%2C%20gene%20correlations%20and%20radar%20charts/RadarCharts_FCgeneExpression.R) is the code to perform the radar charts shown in the manuscript.

Data necessary to perform the analyses should be downloaded from Figshare. Please, download and uncompress the [files](https://figshare.com/s/ac59e7686a48faa43938) inside a folder named **_data_**.


## 3. Gene set enrichment analysis

The scripts for this analysis can be found in this [repository](https://github.com/miquelduranfrigola/sysmalvac_enrichment).


## 4. Cell phenotyping

   * [`Phenotyping_analysis.R`](https://github.com/gmoncu/SysMalVac/blob/master/4_Cell%20phenotyping/Phenotyping_analysis.R) performs the analysis and figures included in the mansucript.
- **_data_** folder: 
   * [`pheno_all.csv`](https://github.com/gmoncu/SysMalVac/blob/master/4_Cell%20phenotyping/data/pheno_all.csv) is a database with the flow cytometry phenotyping and clinical study variables.
