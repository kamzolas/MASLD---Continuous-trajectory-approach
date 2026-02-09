library(igraph) # 2.2.1
library(jsonlite) # 2.0.0
source('../library.R')


################################################################################
# Description
################################################################################
# Enrichment pathway (Reactome) analysis on the consensus network sigantures.
# Output files:
# - A tsv file for each sliding window and direction of differentiation.
################################################################################


################################################################################
# Inputs
################################################################################
main_dir = paste('../../results/ucam_sanyal/networks/', sep='')
load(paste(main_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
load(paste(main_dir, 'network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight.RData', sep=''))
load('../../data/annotation_databases/annotation_lists.RData')
output_dir <-paste(main_dir, 'network_analysis/Reactome_for_consensus_propagation_signatures/', sep='')


################################################################################
# Enrichment analysis
################################################################################
network <- MASLD_unified_undirected_network
ref_annotation_list <- annotation_lists[["Reactome_2024"]]
background <- unique(unname(unlist(ref_annotation_list)))
dataset_background <- as_ids(V(network))
gene_set_for_adjustment <- intersect(dataset_background, background)
length(gene_set_for_adjustment)

dir.create(output_dir, showWarnings = FALSE)
for (sw in names(propagation_signatures)) {
  for (direction in names(propagation_signatures[[sw]])) {
    genes_set <- propagation_signatures[[sw]][[direction]]$genes_signature
    res = execute_fisher_exact_test(genes_set, ref_annotation_list, TRUE, gene_set_for_adjustment)
    res <- res[order(res$adj_p.value, decreasing = FALSE),]
    filename <- paste(output_dir, sw, '_', direction, '.tsv', sep='')
    write.table(res, file = filename, quote=FALSE, sep='\t')
  }
}


