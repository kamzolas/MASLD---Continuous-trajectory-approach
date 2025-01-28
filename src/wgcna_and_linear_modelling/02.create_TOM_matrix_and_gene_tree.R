suppressMessages(library(WGCNA)) # 1.72.1 (enableWGCNAThreads, adjacency, TOMsimilarity)
suppressMessages(library(flashClust)) # 1.1.2 (flashClust)


output_dir <- '../../results/wgcna_and_linear_modelling/'

################################################################################
# Set up the parallel calculation back-end
################################################################################
enableWGCNAThreads()


################################################################################
# Load the objects from the previous step
################################################################################
load(paste(output_dir, "soft_threshold_results.RData", sep=''))


################################################################################
# Calculate the network adjacency and TOM distance matrices
################################################################################
softPower = 4 # found as the optimal value in the previous step
adjacency_matrix = adjacency(GeneXData, power = softPower, type = "signed")
TOM = TOMsimilarity(adjacency_matrix, TOMType="signed")
TOM_distance_matrix = 1-TOM
gene_tree = flashClust(as.dist(TOM_distance_matrix), method="average")


################################################################################
# Save the TOM distance matrix, as well as the hierarchical tree
################################################################################
save(TOM_distance_matrix, file=paste(output_dir, "TOM_distance_matrix.RData", sep=''))
save(gene_tree, file=paste(output_dir, "gene_tree.RData", sep=''))

