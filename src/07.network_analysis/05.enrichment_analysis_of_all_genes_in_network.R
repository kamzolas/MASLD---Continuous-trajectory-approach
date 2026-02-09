library(igraph) # 2.2.1
library(jsonlite) # 2.0.0
source('../library.R')

################################################################################
# Description
################################################################################
# Enrichment pathway (Reactome) analysis on the whole MASLD network
# Output files:
# - A tsv file with the results of pathway analysis
################################################################################


################################################################################
# Inputs
################################################################################
main_dir = paste('../../results/ucam_sanyal/networks/', sep='')
load(paste(main_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
load('../../data/annotation_databases/annotation_lists.RData')


################################################################################
# Enrichment analysis
################################################################################
network <- MASLD_unified_undirected_network
genes_set <- as_ids(V(network))
ref_annotation_list <- annotation_lists[["Reactome_2024"]]
background <- unique(unname(unlist(ref_annotation_list)))
res = execute_fisher_exact_test(genes_set, ref_annotation_list, FALSE)
res <- res[order(res$adj_p.value, decreasing = FALSE),]
filename <- paste(main_dir, 'MASLD_nodes_Reactome.tsv', sep='')
write.table(res, file = filename, quote=FALSE, sep='\t')


