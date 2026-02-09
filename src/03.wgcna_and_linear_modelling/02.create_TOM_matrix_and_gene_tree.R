suppressMessages(library(WGCNA)) # 1.73.0 (enableWGCNAThreads, adjacency, TOMsimilarity)
suppressMessages(library(flashClust)) # 1.1.2 (flashClust)

################################################################################ 
# Description
################################################################################
# Construction of the topological overlap matrix (TOM) and the hierarchical
# clustering of genes given the soft_threshold_results.RData file which contains
# all the necessary data, derived from the previous step. The user should define 
# as input the optimal soft-threshold power.
# Outputs:
# - TOM_distance_matrix.RData
# - gene_tree.RData
################################################################################


################################################################################ 
# Inputs
################################################################################
# Main directory where the results from soft-threshold detection have been saved.
# In this directory two RData object will be saved after the execution of this
# script.
main_dir <- '../../results/ucam_sanyal/wgcna_and_linear_modelling/'
# The optimal value of soft power detected in the previous step 
softPower = 4 # found from the analysis in 01.wgcna_soft_threshold_detection.R


################################################################################
# 1. Set up the parallel calculation back-end
################################################################################
enableWGCNAThreads()


################################################################################
# 2. Load the objects from the previous step
################################################################################
load(paste(main_dir, "soft_threshold_results.RData", sep=''))


################################################################################
# 3. Calculate the network adjacency and TOM distance matrices
################################################################################
adjacency_matrix = adjacency(GeneXData, power = softPower, type = "signed")
TOM = TOMsimilarity(adjacency_matrix, TOMType="signed")
TOM_distance_matrix = 1-TOM
gene_tree = flashClust(as.dist(TOM_distance_matrix), method="average")


################################################################################
# 4. Save the TOM distance matrix, as well as the hierarchical tree
################################################################################
save(TOM_distance_matrix, file=paste(main_dir, "TOM_distance_matrix.RData", sep=''))
save(gene_tree, file=paste(main_dir, "gene_tree.RData", sep=''))

