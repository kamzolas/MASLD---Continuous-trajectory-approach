# MASLD---Continuous-trajectory-approach
By moving beyond conventional stage-based MASLD classification, we uncover the sequence of critical molecular events that enhance our understanding of MASLD/MASH pathophysiology. This approach enabled the identification of novel trajectory-specific biomarkers, offering a more refined and personalised strategy for managing MASLD patients.


## Code archtiecture

The different parts of the analysis have been organised in specific sub-folders in the main 'src' folder.

1. **finding_of_optimal_SWs_sequence**: contains the scripts used to find the optimal set of Sliding Windows for patients' stratification in the pseudo-temporal space.
2. **wgcna_and_linear_modelling**: contains the scripts used to identify the gene co-expression modules (running WGCNA) and to model the MASLD phenotypes using the calcualted eigen-genes.
