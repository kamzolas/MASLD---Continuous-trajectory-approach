# MASLD---Continuous-trajectory-approach
By moving beyond conventional stage-based MASLD classification, we uncover the sequence of critical molecular events that enhance our understanding of MASLD/MASH pathophysiology. This approach enabled the identification of novel trajectory-specific biomarkers, offering a more refined and personalised strategy for managing MASLD patients.


## Code structure

The different parts of the analysis have been organised in specific files and/or sub-folders in the main 'src' folder.

1. File organisation, preprocessing, and main trajectory analysis.
   a. Normalisation and PCA,
   b. Batch correction,
   c. Trajectory analysis and visualisation
2. **finding_of_optimal_SWs_sequence**: contains the scripts used to find the optimal set of Sliding Windows for patients' stratification in the pseudo-temporal space.
3. **deseq_and_msviper**: contains the scripts that execute the differential expression and transcription factor analysis for both discrete and SW-based patients' stratification.
4. **wgcna_and_linear_modelling**: contains the scripts used to identify the gene co-expression modules (running WGCNA) and to model the MASLD phenotypes using the calcualted eigen-genes.
5. **enrichment_analysis**: contains the functions used to perform TF- and pathway (Reactome) enrichment analysis in order to connect the gene co-expression modules with TFs and upstream signalling pathways.
6. **create_MASLD_network**: contains all the steps used to create the MASLD network which constitute the main framework used for the downstream analysis. Firstly, a reference network is created based on available signalling and metabolic pathway databases. Then the individual upstream and downstream networks that are associated with the MASLD variables are constructed. Finally, these distinct networks are unified into the global MASLD network, whose edges are annotated with the functional similarities of interacting pairs.