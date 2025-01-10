library(igraph)
library(rlang)
library(jsonlite)

remove_deadend_nodes <- function(sw_graph) {
  # remove nodes without incoming edges
  in_degrees <- degree(sw_graph, mode='in')
  to_remove_nodes <- names(in_degrees[in_degrees == 0])
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'start']
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'end']
  sw_graph = delete.vertices(sw_graph, v=to_remove_nodes)
  # remove nodes without outcoming edges
  out_degrees <- degree(sw_graph, mode='out')
  to_remove_nodes <- names(out_degrees[out_degrees == 0])
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'end']
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'start']
  sw_graph = delete.vertices(sw_graph, v=to_remove_nodes)
  return(sw_graph)
} 



filter_similar_incoming_edges <- function(sw_graph) {
  
  nodes <- as_ids(V(sw_graph))
  nodes <- as_ids(V(sw_graph))[2:(length(nodes)-1)]
  edges_df <- igraph::as_data_frame(sw_graph)
  
  ### remove from the list of nodes those of start and end
  to_remove <- c()
  for (node in nodes) {
    # get the incoming neighbors
    incoming_neighbors <- as_ids(neighbors(sw_graph, v=node, mode='in'))
    if (length(incoming_neighbors) == 1) {
      next()
    }
    # group them according to their initial sample
    groups <- list()
    for (neighbor in incoming_neighbors) {
      tmp_id <- strsplit(neighbor, split = '_')[[1]][2]
      first_sample <- strsplit(tmp_id, split = '-')[[1]][1]
      groups[[first_sample]] <- c(groups[[first_sample]], neighbor)
    }
    # for each group keep the node with highest deg value or the biggest size
    # remove the incoming edges for all the other nodes
    for (first_sample in names(groups)) {
      neighbors <- groups[[first_sample]]
      if (length(neighbors) == 1) {
        next()
      } else {
        tmp_edges <- edges_df[(edges_df$from %in% neighbors) & (edges_df$to == node),]
        tmp_edges$length <- sapply(tmp_edges$from, function(c) {
          length(strsplit(c, split='-')[[1]])
        })
        tmp_edges <- tmp_edges[order(-tmp_edges$degs, -tmp_edges$length),]
        to_remove <- c(to_remove, rownames(tmp_edges)[2:dim(tmp_edges)[1]])
      }
    }
  }
  ### make the reduction
  to_remove <- as.numeric(to_remove)
  sw_graph <- delete_edges(sw_graph, edges = to_remove)
  sw_graph <- remove_deadend_nodes(sw_graph)
  
  return(sw_graph)
}




remove_edges_with_low_signal <- function(sw_graph) {
  
  tmp_sw_graph = duplicate(sw_graph)
  small_degs <- get.edge.attribute(tmp_sw_graph, 'degs')
  small_degs <- small_degs[!is.na(small_degs)]
  small_degs <- small_degs[small_degs < 100] # 100 is the lower limit
  thr <- quantile(small_degs, 0.5, na.rm=TRUE)
  thr <- as.numeric(thr)
  
  central_edges <- c()
  
  while(TRUE) {
    
    edge_scores_df <- igraph::as_data_frame(sw_graph)
    edge_scores_df$edge_id <- paste(edge_scores_df$from, '|', edge_scores_df$to)
    edge_scores_df$betweenness <- edge_betweenness(sw_graph)
    edge_scores_df <- edge_scores_df[!is.na(edge_scores_df$degs),]
    edge_scores_df <- edge_scores_df[edge_scores_df$degs < thr,]
    
    tmp_edges <- setdiff(edge_scores_df$edge_id, central_edges)
    edge_scores_df <- edge_scores_df[edge_scores_df$edge_id %in% tmp_edges,]
    edge_scores_df <- edge_scores_df[order(edge_scores_df$betweenness),]
    to_remove_edges <- as.numeric(rownames(edge_scores_df)[1])
    
    if (length(rownames(edge_scores_df)) == 0) {
      break
    }  
    tmp_sw_graph <- duplicate(sw_graph)
    tmp_sw_graph <- delete_edges(tmp_sw_graph, edges=to_remove_edges)
    components <- components(tmp_sw_graph)
    main_component <- which(components$csize == max(components$csize))
    # remove nodes out the main_component
    to_remove_nodes <- names(which(components$membership != main_component))
    to_remove_nodes <- to_remove_nodes[!to_remove_nodes %in% c('start', 'end')]
    tmp_sw_graph = delete.vertices(tmp_sw_graph, v=to_remove_nodes)
    
    tmp_sw_graph <- remove_deadend_nodes(tmp_sw_graph)
    # check again the components
    components <- components(tmp_sw_graph)
    if (components$membership['start'] == components$membership['end']) {
      sw_graph <- duplicate(tmp_sw_graph)
    } else {
      central_edges <- c(central_edges, c(edge_scores_df[1,'edge_id']))
    }
  }
  
  sw_graph <- remove_deadend_nodes(sw_graph)
  return(sw_graph)
  
}


output_dir <- '../../results/finding_of_optimal_SWs_sequence/'
directory <- 'sw_graphs_8_24_03'

load(paste(output_dir, directory, '/', 'sw_graph_final.RData', sep=''))
sw_graph_final <- duplicate(sw_graph)

sw_graph_final <- filter_similar_incoming_edges(sw_graph_final)
sw_graph_final <- remove_edges_with_low_signal(sw_graph_final)
paths <- all_simple_paths(sw_graph_final, from='start', to='end', mode='out')

length(paths)


optimal_path_key <- 1276
path <- paths[[optimal_path_key]]
path_nodes <- as_ids(path)[2:(length(path)-1)] # remove start and end nodes
sw_data <- list()
sw = 1
for (node in path_nodes) {
  ids_str = strsplit(node, split='_')[[1]][2]
  ids = as.numeric(strsplit(ids_str, split='-')[[1]])
  samples <- rownames(ref_pc_sorting_df[ids,])
  sw_data[[as.character(sw)]] <- samples
  sw = sw + 1
}

sw_df <- data.frame()
for (sw in names(sw_data)) {
  tmp_vector <- paste(sw_data[[sw]], collapse=';')
  sw <- paste('SW', sw, sep = '_')
  tmp_vector <- c(sw, tmp_vector)
  sw_df <- rbind(sw_df, tmp_vector)
}
colnames(sw_df) <- c('sw', 'samples')
write.table(sw_df, file=paste(output_dir, 'sw_samples.csv', sep=''), sep=',', 
            row.names = FALSE, quote = FALSE)
write.table(sw_df, file=paste('../../data/sw_samples.csv', sep=''), sep=',', 
            row.names = FALSE, quote = FALSE)

