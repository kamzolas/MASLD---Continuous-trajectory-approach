library(igraph) # 1.4.3
library(tidyr) # 1.3.1
library(dplyr) # 1.1.4



################################################################################
#
# Output files:
#
# Two new RData files with the updates version of downstream and upstream 
# networks (labelled as v1), in the 'MASLD_network_construction' folder:
#
# - MASLD_directed_downstream_networks_v1.RData
# - MASLD_directed_downstream_networks_v1.RData
#
################################################################################



################################################################################
#
# The main worker to integrate the semantic measures in network edges. The three
# generated semantic similarities matrices are combined to create a unique score
# for each interaction. Interaction for which there is not any score in the
# matrices are removed. Interaction which contained CHEBI nodes are scored with 
# 1. Finally this process ends up to create a smaller networks, as some nodes 
# have been removed due to the lack of GO annotation.
#
################################################################################
integrate_semantics_worker <- function(input_network) {

    edges_df <- igraph::as_data_frame(input_network)
    
    from_chebi <- which(grepl(x = edges_df$from, pattern='^CHEBI'))
    to_chebi <- which(grepl(x = edges_df$to, pattern='^CHEBI'))
    chebi_indexes <- union(from_chebi, to_chebi)
    edges_df$edge_id <- row.names(edges_df)
    
    ## Load the semantic similarity matrices
    filenames <- list.files(paste(main_dir, 'semantic_similarities/', sep=''))
    for (filename in filenames) {
      sims = read.csv(paste(main_dir, 'semantic_similarities/', filename, sep=''), 
                      sep=',', header = TRUE, row.names = 1, check.names = FALSE)
      measure <- gsub('.csv', '', filename)
      #print(measure)
      #print(dim(sims))
      edges_df[,measure] <- sapply(1:dim(edges_df)[1], function(index) {
        if (index %in% chebi_indexes) {
          sim = 1
        } else {
          row <- edges_df[index,]
          x <- as.character(row$from)
          y <- as.character(row$to)
          sim <- as.numeric(sims[x, y])
          if (length(sim) == 0) {
            sim <- as.numeric(sims[y, x])
            if (length(sim) == 0) {
              sim = NA
            }
          }
        }
        sim
      }, simplify = TRUE)
    }
    
    # Average the semantic similarity scores for each edge
    edges_df$weight <- rowMeans(edges_df[,c('Lin_BP', 'Resnik_BP', 'Wang_BP')], na.rm = TRUE)
    
    # Generate a new network, identical to the input one, but with an additional
    # edge attribute for the semantics
    tmp_network <- set_edge_attr(input_network, 'weight',
                                 index=as.integer(edges_df$edge_id),
                                 value=as.numeric(edges_df$weight))
    
    # Find the edges with NaN semantic similarities and delete them from the new
    # network
    to_remove <- which(is.nan(edges_df$weight))
    to_remove_edges <- as.integer(edges_df[to_remove,'edge_id'])
    tmp_network <- delete_edges(tmp_network, to_remove_edges)
    
    # Check if the derived network is not connected
    components <- components(tmp_network)
    biggest_cluster_id <- which.max(components$csize)
    keep_nodes <- V(tmp_network)[components$membership == biggest_cluster_id]
    output_network <- igraph::induced_subgraph(tmp_network, keep_nodes)
    
    print(paste('duplicates', sum(duplicated(igraph::as_data_frame(output_network)[,c('from', 'to')]))))
    print(paste(is_directed(input_network), is_directed(output_network)))
    print(is_connected(output_network))
    
    return(output_network)
}



################################################################################
#
# 1. Input arguments
#
################################################################################
main_dir = '../../results/networks/'



################################################################################
#
# 2. Transformation of the downstream networks (v0). Their new versions are 
# labelledas v1.
#
################################################################################
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_downstream_networks_v0.RData', sep='')
load(filename)
tmp_list <- list()
for (name in names(downstream_networks_per_var)) {
  network <- integrate_semantics_worker(downstream_networks_per_var[[name]])
  print(paste(name, length(downstream_networks_per_var[[name]]), length(network)))
  tmp_list[[name]] <- network
}
downstream_networks_per_var <- tmp_list
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_downstream_networks_v1.RData', sep='')
save(downstream_networks_per_var, file=filename)
#[1] "Fibrosis 2261 2136"
#[1] "Steatosis 2303 2162"
#[1] "Ballooning 1853 1756"
#[1] "NAS 2151 2033"
#[1] "diff_exp 1044 1024"



################################################################################
#
# 2. Transformation of the upstream networks (v0). Their new versions are 
# labelled as v1.
#
################################################################################
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_upstream_networks_v0.RData', sep='')
load(filename)
tmp_list <- list()
for (name in names(upstream_networks_per_var)) {
  network <- integrate_semantics_worker(upstream_networks_per_var[[name]])
  print(paste(name, length(upstream_networks_per_var[[name]]), length(network)))
  tmp_list[[name]] <- network
}
upstream_networks_per_var <- tmp_list
filename <- paste(main_dir, 'MASLD_network_construction/MASLD_directed_upstream_networks_v1.RData', sep='')
save(upstream_networks_per_var, file=filename)
#[1] "Fibrosis 729 695"
#[1] "Steatosis 306 295"
#[1] "Ballooning 791 752"
#[1] "NAS 719 683"
#[1] "diff_exp 767 730"


















