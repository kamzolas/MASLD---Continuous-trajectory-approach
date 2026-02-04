# MASLD---Continuous-trajectory-approach
By moving beyond conventional stage-based MASLD classification, we uncover the sequence of critical molecular events that enhance our understanding of MASLD/MASH pathophysiology. This approach enabled the identification of novel trajectory-specific biomarkers, offering a more refined and personalised strategy for managing MASLD patients.


## Code structure

The different parts of the analysis have been organised in specific files and/or sub-folders in the main 'src' folder.

1. **preprocessing_and_trajectory_analysis**: File organisation, preprocessing (normalisation, PCA, batch correction) and main trajectory inference and analysis.
**Validation:** The following subfolders include the scripts used for validation of the preprocessing and trajectory inference part of the analysis. 
**a)Linking_disease_state_to_phenotypes**: Script producing the plots showing agreement of progressing disease characteristics with trajectory position, 
**b)Sex_Differences**: Script producing the correlation plots (disease and general population characteristics against the first 10 principal components) separately for males and females, 
**c)Error_bars_in_histology_per_SW**: Script producing the plots of the average histological scores per SW (including the error bars in each SW), 
**d)Correlation_analysis_normalised_vs_batchcorrected_counts**: Script producing all the supplementary plots validating successful batch effect correction, including correlations per disease stage and stratified correlation (all disease stages)
3. **finding_of_optimal_SWs_sequence**: Contains the scripts used to find the optimal set of Sliding Windows for patients' stratification in the pseudo-temporal space.
4. **deseq_and_msviper**: Contains the scripts that execute the differential expression and transcription factor analysis for both discrete and SW-based patients' stratification.
5. **wgcna_and_linear_modelling**: Contains the scripts used to identify the gene co-expression modules (running WGCNA) and to model the MASLD phenotypes using the calcualted eigen-genes.
6. **enrichment_analysis**: Contains the functions used to perform TF- and pathway (Reactome) enrichment analysis in order to connect the gene co-expression modules with TFs and upstream signalling pathways.
7. **create_MASLD_network**: Contains all the steps used to create the MASLD network which constitute the main framework used for the downstream analysis. Firstly, a reference network is created based on available signalling and metabolic pathway databases. Then the individual upstream and downstream networks that are associated with the MASLD variables are constructed. Finally, these distinct networks are unified into the global MASLD network, whose edges are annotated with the functional similarities of interacting pairs.
8. **network_analysis**: Includes the functions used for network propagation on the MASLD network in order to extract the "up" and "down" signatures for each disease stage, as well as the interpretation of these signatures using the Reactome pathways database.
9. **cell_type_deconvolution**: Contains a simple script to aggregate the results from cell type deconvultion and calculate an average composition for each cell type.
10. **pseudo_bulk_analysis**: This is a worklfow of different tasks to generate pseudo bulk RNA-seq samples for each disease state, given the respective cell type proportions and then perform differential and network analysis to identify molecular changes along the disease trajectory. It is assumed that these changes have been derived only from the differentiations in cell type proportions. Finally the part of differentiated pathways that cannot be explained by this factor is defined as the 'deregulated' disease component and is associated with the different cell types to pinpoint disease-relevant cell type-specific processes.
11. **markers_analysis**: Contains scripts to extend the signature of 57 plasma markers with genes functionally related to them and deregulated in the same disease stage. Additional ChEMBL API is used to retrieve drug-treatment data for these markers.
12. **RandomForest_analysis**: Includes all Random Forest classification and regression analyses, applied to produce the biomarker sets described in the study. The training has been applied exclusively on the UCAM/VCU dataset. Additionally, we provide all the files that reproduce the analyses on the external validation transcriptomic datasets (EPoS, GUBRA, Fujiwara) as well as the validation on the plasma proteomic datasets.
13. **Genetic_evidence_57biomarkers**: Contains the script that produces all files and figures related to the GWAS catalog analysis including the gene/biomarker-trait category associations, and enrichment of MASLD biomarker genes on the trait categories (e.g., metabolic, liver, cardiovascular, inflammatory etc), as well as summary gwas files and associated genes/biomarkers with multiple gwas categories.

* Throughout all those scripts, comments have been used to navigate first-time users to run and reproduce the results of our analysis. 

## Software requirements


#### Python libraries:

Python3: 3.12.3

- chembl_webresource_client: 0.10.9 (Python 3.8 environment)
- igraph: 0.11.4
- json: 2.0.9
- numpy: 1.26.4
- pandas: 2.1.4
- scipy: 1.11.4
- statmodels: 0.14.1



#### R libraries:

R: 4.4.1

- assertr: 3.0.1
- BiocParallel: 1.40.2
- biomaRt: 2.54.0
- broom: 1.0.10
- caret: 6.0-94
- cowplot: 1.1.1
- data.table: 1.17.8
- decoupleR: 2.12.0
- DESeq2: 1.46.0
- dplyr: 1.1.4
- flashClust: 1.1.2
- ggplot2: 4.0.1
- ggpubr: 0.6.0
- GO.db: 3.20.0
- GOSemSim: 2.29.2
- gridExtra: 2.3
- igraph: 2.2.1
- jsonlite: 2.0.0
- modelr: 0.1.11
- numbers: 0.8.5
- OmnipathR: 3.8.0
- PCSF: 0.99.1
- pheatmap: 1.0.12
- plotly: 4.11.0
- preprocessCore: 1.68.0
- pROC: 1.18.5
- purrr: 1.2.0
- randomForest: 4.6.14
- Rcpp: 1.1.0
- readr: 2.1.6
- readxl: 1.4.5
- reshape2: 1.4.4
- rlang: 1.1.6
- Seurat: 5.3.1
- SingleCellExperiment: 1.28.1
- slingshot: 3.22
- stats: 4.0.3
- stringr: 1.6.0
- sva: 3.38.0
- tidyr: 1.3.1
- tidyverse: 2.0.0
- viper: 1.40.0
- WGCNA: 1.73.0
  
