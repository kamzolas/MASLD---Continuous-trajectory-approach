library(igraph)
library(jsonlite)


################################################################################
#
# 1. Input arguments
#
################################################################################
key_dir = '../../results/pseudo_bulk_analysis/networks/'
tmp_dir1 = paste(key_dir, 'network_analysis/propagation_results/', sep='')
tmp_dir2 = paste(key_dir, 'network_analysis/consensus_propagation_signatures/', sep='')

files <- list.files(paste(key_dir, 'network_analysis/propagation_results/', sep=''))
files <- files[grepl('.RData', files)]





################################################################################
#
# 2. MASLD network and propagation results loading
#
################################################################################
# Undirected MASLD network
main_project_dir <- '../../results/networks/'
load(paste(main_project_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- MASLD_unified_undirected_network


################################################################################
#
# 3. Filtering and saving the results from network propagation
#
################################################################################
dir.create(tmp_dir2, showWarnings = FALSE, recursive = TRUE)

propagation_signatures <- list()

for (file in files) {
  tmp <- strsplit(x=file, split = '\\.')[[1]][1]
  param <- paste(strsplit(x=tmp, split = '-')[[1]][1], strsplit(x=tmp, split = '-')[[1]][2], sep='-')
  damping <- strsplit(x=tmp, split = '-')[[1]][3]
  filename <- paste(tmp_dir1, file, sep='')
  load(filename)
  tmp_signatures <- list()
  directions <- c('up', 'down')
  for (sw in names(propagation_results)) {
    for (direction in directions) {
      pvalues_df <- propagation_results[[sw]][[direction]][['pvalues']]
      nodes_df <- pvalues_df[pvalues_df$pvalue < 0.05,,drop=FALSE]
      propagated_nodes <- rownames(nodes_df)
      tmp_signatures[[sw]][[direction]]<- propagated_nodes
    }
  }
  propagation_signatures[[param]][[damping]] <- tmp_signatures
}



final_signatures <- list()
final_signatures_json <- list()

for (param in names(propagation_signatures)) {
  
  tmp_signatures <- list()
  
  for (damping in names(propagation_signatures[[param]])){
    tmp <- propagation_signatures[[param]][[damping]]
    for (sw in names(tmp)) {
      tmp_signatures[[sw]][['up']] <- c(tmp_signatures[[sw]][['up']], tmp[[sw]][['up']])
      tmp_signatures[[sw]][['down']] <- c(tmp_signatures[[sw]][['down']], tmp[[sw]][['down']])
    }
  }
  
  for (sw in names(tmp_signatures)) {
    
    # Unify and filter the up_genes
    up_counts <- table(tmp_signatures[[sw]][['up']])
    up_genes <- names(up_counts[up_counts == 3])
    propagated_network <- induced_subgraph(network, vids=up_genes, impl="create_from_scratch")
    detected_components <- components(propagated_network)
    components <- which(detected_components$csize >=3)
    component_nodes <- names(detected_components$membership[detected_components$membership %in% components])
    propagated_network <- induced_subgraph(propagated_network, vids=component_nodes, 
                                           impl="create_from_scratch")
    final_signatures[[param]][[sw]][['up']][['genes_signature']] <- component_nodes
    final_signatures[[param]][[sw]][['up']][['subnetwork']] <- propagated_network
    final_signatures_json[[param]][[sw]][['up']] <- component_nodes
    
    # Unify and filter the down_genes
    down_counts <- table(tmp_signatures[[sw]][['down']])
    down_genes <- names(down_counts[down_counts == 3])
    propagated_network <- induced_subgraph(network, vids=down_genes, impl="create_from_scratch")
    detected_components <- components(propagated_network)
    components <- which(detected_components$csize >=3)
    component_nodes <- names(detected_components$membership[detected_components$membership %in% components])
    propagated_network <- induced_subgraph(propagated_network, vids=component_nodes, 
                                           impl="create_from_scratch")
    final_signatures[[param]][[sw]][['down']][['genes_signature']] <- component_nodes
    final_signatures[[param]][[sw]][['down']][['subnetwork']] <- propagated_network
    final_signatures_json[[param]][[sw]][['down']] <- component_nodes
  }
} 


for (param in names(final_signatures)) {
  #RData
  propagation_signatures <- final_signatures[[param]]
  filename <- paste(tmp_dir2, param, '.RData', sep='')
  save(propagation_signatures, file=filename)
  # JSON
  propagation_signatures_json <- final_signatures_json[[param]]
  filename <- paste(tmp_dir2, param, '.json', sep='')
  write_json(propagation_signatures_json, filename, auto_unbox = TRUE, pretty = FALSE)
}



