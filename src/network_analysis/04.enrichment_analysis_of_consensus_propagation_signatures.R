library(igraph)
library(jsonlite)
source('../library.R')

main_dir = paste('../../results/networks/', sep='')

load(paste(main_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- MASLD_unified_undirected_network


################################################################################
# Enrichment analysis of network
################################################################################
load('../../data/annotation_databases/annotation_lists.RData')
ref_annotation_list <- annotation_lists[["Reactome_2024"]]
background <- unique(unname(unlist(ref_annotation_list)))
dataset_background <- as_ids(V(network))
gene_set_for_adjustment <- intersect(dataset_background, background)
length(gene_set_for_adjustment)

load(paste(main_dir, 'network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight.RData', sep=''))

output_dir <-paste(main_dir, 'network_analysis/Reactome_for_consensus_propagation_signatures/', sep='')
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






