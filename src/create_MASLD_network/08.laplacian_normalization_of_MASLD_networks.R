library(igraph) # 1.4.3
source("../library.R")


################################################################################
#
# Output files:
#
# The two formats of unified network with normalized edge weights based on 
# laplacian normalization:
#
# - MASLD_unified_undirected_network_with_semantics.RData
# - MASLD_unified_directed_network_with_semantics.RData
#
################################################################################



################################################################################
#
# 1. Input arguments
#
################################################################################
main_dir = '../../results/networks/'



################################################################################
#
# 2. Loading and transformation of the undirected network.
#
################################################################################
load(paste(main_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- laplacian_normalization(MASLD_unified_undirected_network)
filename <- paste(main_dir, "MASLD_unified_undirected_network_with_semantics.RData", sep='')
MASLD_unified_undirected_network <- network
save(MASLD_unified_undirected_network, file=filename)



################################################################################
#
# 2. Loading and transformation of the directed network.
#
################################################################################
load(paste(main_dir, 'MASLD_unified_directed_network_with_semantics.RData', sep=''))
network <- laplacian_normalization(MASLD_unified_directed_network)
filename <- paste(main_dir, "MASLD_unified_directed_network_with_semantics.RData", sep='')
MASLD_unified_directed_network <- network
save(MASLD_unified_directed_network, file=filename)



