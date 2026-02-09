library(igraph) # 2.2.1
library(rlang) # 1.1.6
library(jsonlite) # 2.0.0

################################################################################ 
# Description
################################################################################
# Filtering of the sw_graph object generated in the previous step. The graph is
# filtered in order to reduce the amount of potential SW-based trajectories. 
# The selection of the optimal used is performed in the next step.
# Outputs:
# - A json file the includes all the potential SW-based trajectories (final_paths_after_filtering.json)
# - Igraph objects which contain all the potential SW-based trajectories (final_paths_in_sw_graph.RData)
# - A tsv file which contains some measurements about the potential paths. These values
# will be used to find the optimal paths (final_paths_stats.tsv).
################################################################################


################################################################################
# Function to remove nodes nodes without incoming or outcoming edges from 
# sw_graph
################################################################################
remove_deadend_nodes <- function(sw_graph) {
  # remove nodes without incoming edges
  in_degrees <- degree(sw_graph, mode='in')
  to_remove_nodes <- names(in_degrees[in_degrees == 0])
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'start']
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'end']
  sw_graph = delete_vertices(sw_graph, v=to_remove_nodes)
  # remove nodes without outcoming edges
  out_degrees <- degree(sw_graph, mode='out')
  to_remove_nodes <- names(out_degrees[out_degrees == 0])
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'end']
  to_remove_nodes <- to_remove_nodes[to_remove_nodes != 'start']
  sw_graph = delete_vertices(sw_graph, v=to_remove_nodes)
  return(sw_graph)
} 


################################################################################
# Function to filter sw_graph edges when they lead to the same target node from 
# different sources. For such a group of edges, the edge with highest deg value 
# or that which originates from the source node with the biggest number of 
# samples will remain on the sw_graph.
################################################################################
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


################################################################################
# Function to filter out low signal edges keeping the connectivity of the 
# sw_graph.
################################################################################
remove_edges_with_low_signal <- function(sw_graph) {
  
  tmp_sw_graph = duplicate(sw_graph)
  small_degs <- edge_attr(tmp_sw_graph, 'degs')
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
    tmp_sw_graph = delete_vertices(tmp_sw_graph, v=to_remove_nodes)
    
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


################################################################################
# Main process
################################################################################
output_dir <- '../../results/ucam_sanyal/finding_of_optimal_SWs_sequence/'
sorted_samples_filename <- "../../data/ucam_sanyal/PC1_sorted_samples.csv"
sorted_samples_df = read.table(sorted_samples_filename, sep=',', header=TRUE)
rownames(sorted_samples_df) = sorted_samples_df[,1] 
sorted_samples_df[,1] = NULL
sorted_samples_df[,2] = rownames(sorted_samples_df)
colnames(sorted_samples_df) <- c('Sorting_axis', 'Sample_name')

directories <- c('sw_graph_6_15_02', 'sw_graph_8_17_025', 'sw_graph_8_24_03',
                 'sw_graph_10_20_03', 'sw_graph_12_25_035')

paths_data <- list()
final_graphs <- list()

for (directory in directories) {
  load(paste(output_dir, directory, '/', 'sw_graph_final.RData', sep=''))
  sw_graph_final <- duplicate(sw_graph)

  print(paste('initial number of nodes:', length(V(sw_graph_final))))
  print(paste('initial number of edges:', length(E(sw_graph_final))))

  # Initial paths
  # Be careful!
  #paths <- all_simple_paths(sw_graph_final, from='start', to='end', mode='out')
  #print(paste('initial length of paths:', length(paths)))
  
  ##############################################################################
  # Filtering 1: Concatenate incoming edges 
  ##############################################################################
  sw_graph_final <- filter_similar_incoming_edges(sw_graph_final)
  
  print(paste('filtering1 number of nodes:', length(V(sw_graph_final))))
  print(paste('filtering1 number of edges:', length(E(sw_graph_final))))
  
  paths <- all_simple_paths(sw_graph_final, from='start', to='end', mode='out')
  print(paste('filtering1 length of paths:', length(paths)))
  
  #components <- components(sw_graph_final)
  is_connected(sw_graph_final) # shoud be true
  
  ##############################################################################
  # Filtering 2: Remove edges edges with deg values below the 50th quantile
  ##############################################################################
  sw_graph_final <- remove_edges_with_low_signal(sw_graph_final)
  print(paste('filtering2 number of nodes:', length(V(sw_graph_final))))
  print(paste('filtering2 number of edges:', length(E(sw_graph_final))))
  
  #components <- components(sw_graph_final)
  is_connected(sw_graph_final) # shoud be true
  
  paths <- all_simple_paths(sw_graph_final, from='start', to='end', mode='out')
  print(paste('filtering2 number of paths:', length(paths)))
  
  paths_data[[directory]] <- paths
  final_graphs[[directory]] <- sw_graph_final
  
  print('END')
}


################################################################################
# Get the paths, retrieve the sample names and save the final json file
################################################################################
json_data <- list()
paths_in_sw_graph <- list()

for (directory in names(paths_data)) {
  sw_graph <- final_graphs[[directory]]
  j = 0
  for (path in paths_data[[directory]]) {
    j <- j+1
    path_nodes <- as_ids(path)[2:(length(path)-1)] # remove start and end nodes
    sw_data <- list()
    sw = 1
    for (node in path_nodes) {
      ids_str = strsplit(node, split='_')[[1]][2]
      ids = as.numeric(strsplit(ids_str, split='-')[[1]])
      samples <- rownames(sorted_samples_df[ids,])
      sw_data[[as.character(sw)]] <- samples
      sw = sw + 1
    }
    path_key <- paste(directory, as.character(j), sep='_')
    json_data[[path_key]] <- sw_data
    path_subgraph <- induced_subgraph(sw_graph, vids=path_nodes)
    paths_in_sw_graph[[path_key]] <- path_subgraph
  }
}

filename <- paste(output_dir, 'final_paths_after_filtering.json', sep='')
write_json(json_data, filename)

stats_df <- list()
for (key in names(paths_in_sw_graph)) {
  subgraph <- paths_in_sw_graph[[key]]
  degs <- edge_attr(subgraph, name='degs')
  pass_lower_thr <- sum(degs > 100)
  cv = sd(degs)/mean(degs)
  vector <- c(median(degs), mean(degs), sd(degs), cv, min(degs), max(degs), 
              pass_lower_thr, length(degs)+1)
  stats_df[[key]] <- vector
}

stats_df <- data.frame(stats_df)
stats_df <- t(stats_df)
colnames(stats_df) <- c('median', 'mean', 'sd', 'cv', 
                        'min_value', 'max_value', 'pass_lower_thr', 'N')
save(paths_in_sw_graph, file = paste(output_dir, 'final_paths_in_sw_graph.RData', sep=''))
write.table(stats_df, file = paste(output_dir, 'final_paths_stats.tsv', sep=''), sep = '\t')

