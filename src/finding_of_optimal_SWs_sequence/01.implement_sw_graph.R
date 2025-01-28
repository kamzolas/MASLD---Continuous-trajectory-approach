library(DESeq2)
library(jsonlite)
library(numbers)
library(rlang)
library(igraph)
library(data.table)

library(BiocParallel)
register(MulticoreParam(4))


################################################################################
#
# DESeq worker
#
################################################################################
deseq_worker <- function(sw1_samples, sw2_samples) {
  
  sw1_samples <- setdiff(sw1_samples, sw2_samples)
  
  # Subset the counts matrix based on the input samples
  sw1_cts <- cts[,sw1_samples,drop=FALSE]
  sw2_cts <- cts[,sw2_samples,drop=FALSE]
  
  # Customize the template and cts for sw1_samples
  sw1_template_df <- template_df[sw1_samples, c('Dataset', 'Sex'), drop=FALSE]
  sw1_key <- 'SW1'
  sw1_template_df$sample_id <- paste(sw1_key, seq(1,length(sw1_samples)), sep='_')
  sw1_template_df$condition <- rep(sw1_key, length(sw1_samples))
  colnames(sw1_cts) = sw1_template_df$sample_id
  
  # Customize the template and cts for sw2_samples
  sw2_template_df <- template_df[sw2_samples, c('Dataset', 'Sex'), drop=FALSE]
  sw2_key <- 'SW2'
  sw2_template_df$sample_id <- paste(sw2_key, seq(1,length(sw2_samples)), sep='_')
  sw2_template_df$condition <- rep(sw2_key, length(sw2_samples))
  colnames(sw2_cts) = sw2_template_df$sample_id
  
  # Concatenate the count matrices and templates
  counts_matrix = cbind(sw1_cts, sw2_cts)
  coldata = rbind(sw1_template_df, sw2_template_df)
  coldata$sample_id == colnames(counts_matrix)
  rownames(coldata) <- NULL
  coldata$condition <- as.factor(coldata$condition)
  
  # Factorization of variables
  coldata$batch_id <- paste0(coldata$Dataset, coldata$Sex)
  coldata$batch_id <- as.factor(coldata$batch_id)
  
  # Run DESeq and get the results
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                  colData = coldata,
                                  design = ~ batch_id + condition)
    dds$condition <- relevel(dds$condition, ref = sw1_key)
    dds <- DESeq(dds, parallel=TRUE)
    res <- results(dds)
    res <- as.data.frame(res)
    res <- res[!is.na(res$padj),]
    degs <- dim(res[res$padj < 0.1,])[1]
  },
  error = function(e) {
    degs <<- 0
  })
  
  return(degs)
}





################################################################################
#
# Function to add progressively new nodes in the graph
#
################################################################################
add_new_nodes <- function(sw_graph) {
  
  # Find the leaf nodes
  if (length(V(sw_graph))==1){
    leaf_nodes <- as_ids(V(sw_graph))
  } else {
    out_degrees <- degree(sw_graph, mode='out')
    leaf_nodes <- as_ids(V(sw_graph))[which(out_degrees == 0)]
  }
  
  # Use leaf nodes as starting points to add new nodes
  for (leaf_node in leaf_nodes) {
    leaf_node_attrs <- vertex_attr(sw_graph, index=leaf_node)
    # Get the starting index (which corresponds to a sample in their sorted list) 
    # for the newcomer nodes, from the attributes of the examined (leaf) node
    start_index <- leaf_node_attrs$start_index
    N <- dim(ref_pc_sorting_df)[1] # total number of samples
    if (start_index == N) {
      next()  # we have reach the end of the graph
    }
    for (size in seq_of_sizes) {
      # Based on the size of the new sample set, define the end index
      end_index <- start_index + size - 1
      if (end_index > N) {
        # if the end index surpass N, then set N as end index
        end_index <- N
      } else if ((N-end_index) < 4) {
        # if the remaining samples are less the 4, then increase the size of the
        # new set by including the remaining samples in it
        end_index <- N
      } else {
        # pass
      }
      # get the new samples and define the corresponding starting index
      samples <- rownames(ref_pc_sorting_df[start_index:end_index,])
      indexes <- which(rownames(ref_pc_sorting_df) %in% samples)
      if (end_index == N) {
        index <- N
      } else {
        index <- (max(indexes)+1) - ceiling(length(samples)*overlap)
      }
      node_id <- paste(paste(indexes, collapse='-'), index, sep='_')
      node_id <- paste('n_', node_id, sep='')
      if (node_id %in% as_ids(V(sw_graph))) {
        ##pass
      } else {
        node_attrs <- list('name' = node_id,
                           'start_index' = index)
        sw_graph <- add_vertices(sw_graph, nv=1, attr = node_attrs)
      }
      sw_graph <- sw_graph + edge(leaf_node, node_id)
    }
  }
  
  return(sw_graph)
}





################################################################################
#
# Node it to samples
#
################################################################################
transform_node_id_to_samples <- function(node_id) {
  samples_str <- strsplit(x = node_id, split='_')[[1]][2]
  samples_id <- strsplit(x = samples_str, split='-')[[1]]
  samples_id <- as.numeric(samples_id)
  samples <- rownames(ref_pc_sorting_df)[samples_id]
  return(samples)
}



################################################################################
#
# Count data and metadata loading
#
################################################################################

ref_pc_sorting_df = read.table("../../data/PC1_sorted_samples.csv", sep=',', header=TRUE)
rownames(ref_pc_sorting_df) = ref_pc_sorting_df[,1] 
ref_pc_sorting_df[,1] = NULL
ref_pc_sorting_df[,2] = rownames(ref_pc_sorting_df)
colnames(ref_pc_sorting_df) <- c('PC1_value', 'Sample.name')

template_df = read.csv(file = "../../data/metadata.csv")[, -c(1,3)]
rownames(template_df) = template_df$Sample.name

cts_ucamsanyal <- read.csv(file = "../../data/merged_counts.csv", sep = ",")
rownames(cts_ucamsanyal) = as.character(cts_ucamsanyal$X)
cts_ucamsanyal = cts_ucamsanyal[,-c(1)]
colnames(cts_ucamsanyal) = template_df$Sample.name
cts_ucamsanyal <- cts_ucamsanyal[(rowSums(cts_ucamsanyal)>dim(cts_ucamsanyal)[2]),] #Exclude low expressed counts
unique(colnames(cts_ucamsanyal) == template_df$Sample.name)
cts = data.matrix(cts_ucamsanyal)

colnames(cts) == row.names(template_df)


################################################################################
#
# Main process
#
################################################################################
#ref_pc_sorting_df <- ref_pc_sorting_df[1:25,]
ref_upper_degs_limit <- 1000
ref_lower_degs_limit <- 100

# 1
#seq_of_sizes <- seq(6,15)
#overlap = 0.2

# 2 
#seq_of_sizes <- seq(8,17)
#overlap = 0.25

# 3
#seq_of_sizes <- seq(10,20)
#overlap <- 0.3

# 4
#seq_of_sizes <- seq(12,25)
#overlap <- 0.35

seq_of_sizes <- seq(15,20)
overlap <- 0.5


################################################################################
# 1. Initialize and build the graph
################################################################################
sw_graph <- make_empty_graph(n=0, directed=TRUE)

# Add the root (start) node
node_attrs <- list('name' = 'start', 'start_index' = 1)
# Add the nodes of the first level and their edges with root
sw_graph <- add_vertices(sw_graph, nv=1, attr = node_attrs)
sw_graph <- add_new_nodes(sw_graph)
current_edges <- length(E(sw_graph))
new_edges <- current_edges

# Add progressively new nodes and edges - built the whole graph
added_edges <- c()
while (new_edges > 0) {
  sw_graph <- add_new_nodes(sw_graph)
  new_edges <- length(E(sw_graph)) - current_edges
  current_edges <- length(E(sw_graph))
  print(paste(length(E(sw_graph)), new_edges))
  added_edges <- c(added_edges, new_edges)
}

# Find the leaves and connect them with the end node 
leaf_nodes <- as_ids(V(sw_graph))[which(degree(sw_graph, mode='out') == 0)]
node_id <- 'end'
node_attrs <- list('name' = 'end', 'start_index' = -1)
sw_graph <- add_vertices(sw_graph, nv=1, attr = node_attrs)
for (leaf_node in leaf_nodes) {
  sw_graph <- sw_graph + edge(leaf_node, node_id)
}

# Remove duplicated edges
sw_graph <- simplify(sw_graph)
print(length(E(sw_graph)))
print(length(V(sw_graph)))
save(sw_graph, file="../../results/finding_of_optimal_SWs_sequence/sw_graph_initial.RData")



################################################################################
# 2. Scan the sw_graph
################################################################################
step = 1
from_nodes <- neighborhood(sw_graph, nodes='start', mode='out', order=step, mindist=step)
from_nodes <- as_ids(from_nodes[[1]])
from_nodes <- from_nodes[from_nodes != 'end']

runs <- c()
while (length(from_nodes) > 0) {
  j <- 0
  examined_edge_ids <- c()
  print(paste('FROM NODES SIZE:', length(from_nodes)))
  degs_distribution <- c()
  
  # Run DESeq for the new edges
  for (from_node in from_nodes) {
    # Retrieve the samples
    sw1_samples <- transform_node_id_to_samples(from_node)
    # Get it's one-step node
    to_nodes <- neighbors(sw_graph, v=from_node, mode='out')$name
    to_nodes <- to_nodes[to_nodes != 'end']
    # For each direct neighbor, get the respective samples and run DESeq
    for (to_node in to_nodes) {
      edge_id <- get.edge.ids(sw_graph, c(from_node, to_node))
      examined_edge_ids <- c(examined_edge_ids, edge_id)
      degs <- edge_attr(sw_graph, name='degs', index=edge_id)
      if (is.na(degs) || is.null(degs)) {
        j <- j + 1
        sw2_samples <- transform_node_id_to_samples(to_node)
        degs <- deseq_worker(sw1_samples, sw2_samples)
        if (is.na(degs)) {
          degs <- 0
        }
        degs_distribution <- c(degs_distribution, c(degs))
        print(paste(from_node, to_node, degs))
        sw_graph <- set_edge_attr(sw_graph, 'degs', edge_id, degs)
      } else {
        degs_distribution <- c(degs_distribution, c(degs))
        #pass
      }
    }
  }
  
  runs <- c(runs, c(j))
  
  degs_distribution <- c(degs_distribution, c(0)) # avoid to have only NAs
  
  if (is.null(degs_distribution)) {
    q_low = ref_lower_degs_limit
    q_up = ref_upper_degs_limit
  } else {
    q_low = quantile(degs_distribution, 0.8, na.rm = TRUE)
    q_up = quantile(degs_distribution, 0.2, na.rm = TRUE)
  }
  if (q_low < ref_lower_degs_limit) {
    lower_degs_limit <- q_low - 1
  } else {
    lower_degs_limit <- ref_lower_degs_limit
  }
  if (q_up > ref_upper_degs_limit) {
    upper_degs_limit <- q_up + 1
  } else {
    upper_degs_limit <- ref_upper_degs_limit
  }
  
  # Remove the edges which do not satisfy the thresholds
  if (is.null(examined_edge_ids)) {
    #pass
  } else {
    examined_edge_ids <- examined_edge_ids[order(examined_edge_ids, decreasing = FALSE)]
    degs_values <- get.edge.attribute(sw_graph, name='degs', index=examined_edge_ids)
    to_remove <- which((degs_values < lower_degs_limit) | (degs_values > upper_degs_limit))
    to_remove_edge_ids <- examined_edge_ids[to_remove]
    sw_graph <- delete_edges(sw_graph, edges = to_remove_edge_ids)
    
    # Remove nodes without connections
    to_remove_nodes <- c('none')
    while (length(to_remove_nodes) > 0) {
      in_degrees <- degree(sw_graph, mode='in')
      to_remove_in_nodes <- names(which(in_degrees == 0))
      to_remove_in_nodes <- to_remove_in_nodes[to_remove_in_nodes != 'start']
      
      out_degrees <- degree(sw_graph, mode='out')
      to_remove_out_nodes <- names(which(out_degrees == 0))
      to_remove_out_nodes <- to_remove_out_nodes[to_remove_out_nodes != 'end']
      to_remove_out_nodes <- to_remove_out_nodes[to_remove_out_nodes != 'start']
      
      to_remove_nodes <- union(to_remove_in_nodes, to_remove_out_nodes)
      sw_graph <- delete_vertices(sw_graph, v = to_remove_nodes)
    }
  }
  
  step = step + 1
  from_nodes <- neighborhood(sw_graph, nodes='start', mode='out',
                             order=step, mindist=step)
  from_nodes <- as_ids(from_nodes[[1]])
  from_nodes <- from_nodes[from_nodes != 'end']
  
  print(paste('GRAPH SIZE:', length(V(sw_graph))))
  filename = paste("sw_graph_tmp_", as.character(step-1), ".RData", sep='')
  save(sw_graph, runs, file=filename)
  
}

save(sw_graph, runs, file="../../results/finding_of_optimal_SWs_sequence/sw_graph_final.RData")
