# MASLD---Continuous-trajectory-approach
By moving beyond conventional stage-based MASLD classification, we uncover the sequence of critical molecular events that enhance our understanding of MASLD/MASH pathophysiology. This approach enabled the identification of novel trajectory-specific biomarkers, offering a more refined and personalised strategy for managing MASLD patients.


## Code structure

The different parts of the analysis have been organised in specific files and/or sub-folders in the main 'src' folder.

1. **preprocessing_and_trajectory_analysis**: File organisation, preprocessing (normalisation, PCA, batch correction) and main trajectory analysis.
2. **finding_of_optimal_SWs_sequence**: contains the scripts used to find the optimal set of Sliding Windows for patients' stratification in the pseudo-temporal space.
3. **deseq_and_msviper**: contains the scripts that execute the differential expression and transcription factor analysis for both discrete and SW-based patients' stratification.
4. **wgcna_and_linear_modelling**: contains the scripts used to identify the gene co-expression modules (running WGCNA) and to model the MASLD phenotypes using the calcualted eigen-genes.
5. **enrichment_analysis**: contains the functions used to perform TF- and pathway (Reactome) enrichment analysis in order to connect the gene co-expression modules with TFs and upstream signalling pathways.
6. **create_MASLD_network**: contains all the steps used to create the MASLD network which constitute the main framework used for the downstream analysis. Firstly, a reference network is created based on available signalling and metabolic pathway databases. Then the individual upstream and downstream networks that are associated with the MASLD variables are constructed. Finally, these distinct networks are unified into the global MASLD network, whose edges are annotated with the functional similarities of interacting pairs.
7. **network_analysis**: includes the functions used for network propagation on the MASLD network in order to extract the "up" and "down" signatures for each disease stage, as well as the interpretation of these signatures using the Reactome pathways database.
8. **cell_type_deconvolution**: contains a simple script to aggregate the results from cell type deconvultion and calculate an average composition for each cell type.
9. **pseudo_bulk_analysis**: this is a worklfow of different tasks to generate pseudo bulk RNA-seq samples for each disease state, given the respective cell type proportions and then perform differential and network analysis to identify molecular changes along the disease trajectory. It is assumed that these changes have been derived only from the differentiations in cell type proportions. Finally the part of differentiated pathways that cannot be explained by this factor is defined as the 'deregulated' disease component and is associated with the different cell types to pinpoint disease-relevant cell type-specific processes.
10. **markers_analysis**: contains scripts to extend the signature of 24 plasma markers with genes functionally related to them and deregulated in the same disease stage. Additional ChEMBL API is used to retrieve drug-treatment data for these markers.


## Software requirements:

#### Python libraries:

Python3: 3.8.10

- pandas: 1.3.5
- numpy: 1.22.4
- json: 2.0.9
- igraph: 0.11.4
- scipy: 1.10.1
- statmodels: 0.14.0
- chembl_webresource_client: 0.10.9

#### R libraries:

R: 4.4.2

- ggplot2: 3.4.2
- plotly: 4.10.4
- tidyverse: 2.0.0
- tidyr: 1.3.0
- dplyr: 1.1.2
- stringr: 1.5.0
- readr: 2.1.4
- readxl: 1.4.3
- jsonlite: 1.8.8
- rlang: 1.1.3
- data.table: 1.14.18
- BiocParallel: 1.37.1
- preprocessCore: 1.65.0
- WGCNA: 1.72.1
- flashClust: 1.1.2
- broom: 1.0.4
- modelr: 0.1.11
- sva: 3.52.0
- slingshot: 2.12.0
- SingleCellExperiment: 1.25.1
- Seurat: 5.2.1
- DESeq2: 1.43.5
- viper: 1.34.0
- purrr: 1.0.1
- numbers: 0.8.5
- igraph: 1.4.3
- OmnipathR: 3.8.0
- PCSF: 0.99.1
- GOSemSim: 2.29.2
- Rcpp: 1.0.12