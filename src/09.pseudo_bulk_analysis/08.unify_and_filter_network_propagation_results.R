library(igraph) # 2.2.1
library(jsonlite) # 2.0.0

################################################################################
# Description
################################################################################
# Network prpagatio nwas executed for different damping factors. In this step, we
# integrate all the results to create a consesus signature for each sliding window
# and differentiation direction. Additioanl filtering is performed to removed 
# isolated nodes from the highly-propagated networks.
# Output files:
# - The consensus sigantures in RData, json and tsv formats
################################################################################


################################################################################
# Inputs
################################################################################
network_dir = '../../results/ucam_sanyal/networks/'
main_dir = '../../results/ucam_sanyal/pseudo_bulk_analysis/networks/'
tmp_dir1 = paste(main_dir, 'network_analysis/propagation_results/', sep='')
tmp_dir2 = paste(main_dir, 'network_analysis/consensus_propagation_signatures/', sep='')

files <- list.files(paste(main_dir, 'network_analysis/propagation_results/', sep=''))
files <- files[grepl('.RData', files)]


################################################################################
# 1. Loading of MASLD network and propagation results
################################################################################
# Undirected MASLD network
load(paste(network_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- MASLD_unified_undirected_network


################################################################################
# 2. Filtering and saving the results from network propagation
################################################################################
dir.create(tmp_dir2, showWarnings = FALSE, recursive = TRUE)

propagation_signatures <- list()
propagation_pvalues <- list()
for (file in files) {
  tmp <- strsplit(x=file, split = '\\.')[[1]][1]
  param <- paste(strsplit(x=tmp, split = '-')[[1]][1], strsplit(x=tmp, split = '-')[[1]][2], sep='-')
  damping <- strsplit(x=tmp, split = '-')[[1]][3]
  filename <- paste(tmp_dir1, file, sep='')
  load(filename)
  tmp_signatures <- list()
  tmp_pvalues <- list()
  directions <- c('up', 'down')
  for (sw in names(propagation_results)) {
    for (direction in directions) {
      pvalues_df <- propagation_results[[sw]][[direction]][['pvalues']]
      tmp_pvalues[[sw]][[direction]] <- pvalues_df
      nodes_df <- pvalues_df[pvalues_df$pvalue < 0.05,,drop=FALSE]
      propagated_nodes <- rownames(nodes_df)
      tmp_signatures[[sw]][[direction]]<- propagated_nodes
    }
  }
  propagation_signatures[[param]][[damping]] <- tmp_signatures
  propagation_pvalues[[param]][[damping]] <- tmp_pvalues
}

final_signatures <- list()
final_signatures_json <- list()

for (param in names(propagation_signatures)) {
  param <- "random_sw-lapl_norm_weight"
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
    components <- which(detected_components$csize >=2)
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
    components <- which(detected_components$csize >=2)
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

for (param in names(propagation_pvalues)) {
  
  directions = c('up', 'down')
  direction_pvalues_list <- list()
  
  for (direction in directions) {
    for (sw in names(propagation_signatures)) {
      tmp_list <- list(propagation_pvalues[[param]][['085']][[sw]][[direction]],
                       propagation_pvalues[[param]][['07']][[sw]][[direction]],
                       propagation_pvalues[[param]][['05']][[sw]][[direction]])
      df <- Reduce(function(x,y) 
        transform(merge(x, y, by.x=0, by.y=0, all=TRUE), row.names='Row.names'), 
        tmp_list)
      df <- -log10(df+0.001)
      df[df < 0] = 0
      df <- as.data.frame(rowMeans(df))
      colnames(df) <- c(sw)
      direction_pvalues_list[[sw]] <- df
    }
    df <- Reduce(function(x,y) 
      transform(merge(x, y, by.x=0, by.y=0, all=TRUE), row.names='Row.names'), 
      direction_pvalues_list)
    filename <- paste(tmp_dir2, param, '_', direction, '.tsv', sep='')
    write.table(df, filename, sep='\t')
  }
}
