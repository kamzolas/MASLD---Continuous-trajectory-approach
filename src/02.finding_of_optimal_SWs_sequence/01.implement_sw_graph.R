library(DESeq2) # 1.46.0
library(jsonlite) # 2.0.0
library(numbers) # 0.8.5
library(rlang) # 1.1.6
library(igraph) # 2.2.1
library(data.table) # 1.17.8
library(assertr) # 3.0.1
library(BiocParallel) # 1.40.2
register(MulticoreParam(4))

################################################################################ 
# Description
################################################################################
# This script is used as the first step to define patient stratification into 
# sliding windows (SWs) (a sequence of overlapping groups [SW1, SW2,..., SWX]). 
# In order to define an optimal SW configuration, a graph-based optimisation 
# framework that evaluates all possible sequences of SWs is implemented here. 
# In this framework, disease stages are represented as nodes in a directed acyclic 
# graph (DAG), where each node corresponds to a group of samples aligned by 
# pseudo-temporal order and edges connect overlapping windows, ensuring continuity 
# along the disease trajectory. The input parameters for this graph-based method 
# is a set of window sizes S=[S1,…, Sm] and a predefined overlap percentage α 
# between adjacent windows. Initially this graph is constructed. Then, for each 
# edge (connection betwwen two SWs) DESeq2 is applied to identify the amount of 
# differentially expressed genes between the two potential disease stages. 
# These amounts of DEGs will be used as criterion to filter the graph, keep 
# the most informative SW paths and select the optimal one in the following 
# scipts.
# Output:
# - The initial graph before DESeq2 implementation sw_graph_initial.RData.
# - The final graph after the implementation of DESeq2 sw_graph_final.RData and
# some filtering.
# - Multiple intermediate versions of the graph (tmp files) just to keep record 
# of the analysis.
################################################################################


################################################################################
# DESeq worker
################################################################################
deseq_worker <- function(sw1_samples, sw2_samples, batches) {
  
  sw1_samples <- setdiff(sw1_samples, sw2_samples)
  
  # Subset the counts matrix based on the input samples
  sw1_counts_matrix <- counts_matrix[,sw1_samples,drop=FALSE]
  sw2_counts_matrix <- counts_matrix[,sw2_samples,drop=FALSE]
  
  # Customize the template and counts_matrix for sw1_samples
  sw1_template_df <- template_df[sw1_samples, batches, drop=FALSE]
  sw1_key <- 'SW1'
  sw1_template_df$sample_id <- paste(sw1_key, seq(1,length(sw1_samples)), sep='_')
  sw1_template_df$condition <- rep(sw1_key, length(sw1_samples))
  colnames(sw1_counts_matrix) = sw1_template_df$sample_id
  
  # Customize the template and counts_matrix for sw2_samples
  sw2_template_df <- template_df[sw2_samples, batches, drop=FALSE]
  sw2_key <- 'SW2'
  sw2_template_df$sample_id <- paste(sw2_key, seq(1,length(sw2_samples)), sep='_')
  sw2_template_df$condition <- rep(sw2_key, length(sw2_samples))
  colnames(sw2_counts_matrix) = sw2_template_df$sample_id
  
  # Concatenate the count matrices and templates
  sw_counts_matrix = cbind(sw1_counts_matrix, sw2_counts_matrix)
  coldata = rbind(sw1_template_df, sw2_template_df)
  coldata$sample_id == colnames(counts_matrix)
  rownames(coldata) <- NULL
  coldata$condition <- as.factor(coldata$condition)
  
  # Factorization of variables
  coldata$batch_id <- as.character(col_concat(coldata[,batches]))
  coldata$batch_id <- as.factor(coldata$batch_id)
  
  # Run DESeq and get the results
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = sw_counts_matrix,
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
# Function to add progressively new nodes in the graph
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
    N <- dim(sorted_samples_df)[1] # total number of samples
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
      samples <- rownames(sorted_samples_df[start_index:end_index,])
      indexes <- which(rownames(sorted_samples_df) %in% samples)
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
# Transform node id to samples names from the dataset
################################################################################
transform_node_id_to_samples <- function(node_id) {
  samples_str <- strsplit(x = node_id, split='_')[[1]][2]
  samples_id <- strsplit(x = samples_str, split='-')[[1]]
  samples_id <- as.numeric(samples_id)
  samples <- rownames(sorted_samples_df)[samples_id]
  return(samples)
}



################################################################################
# Inputs
################################################################################

# 1. A csv file with sample names sorted by PC1 (or any other criterion)
sorted_samples_filename <- "../../data/ucam_sanyal/PC1_sorted_samples.csv"
sorted_samples_df = read.table(sorted_samples_filename, sep=',', header=TRUE)
# Should be like:
#               X          pc1         pc2
#1   GSM3998327  -27.58061891 -0.76788978
#2   GSM3998301  -20.61404364  1.22527002
#3   GSM3998258  -18.34908949 -0.98698704
#4   GSM3998188  -18.27205652  0.91769042
#5   GSM3998332  -17.39840662 -1.22990361
#
rownames(sorted_samples_df) = sorted_samples_df[,1] 
sorted_samples_df[,1] = NULL
sorted_samples_df[,2] = rownames(sorted_samples_df)
colnames(sorted_samples_df) <- c('Sorting_axis', 'Sample_name')

# 2. Samples metadata, the column of sample names and the important batches that
# should be taken into account in differential expression analysis.
template_filename <- "../../data/ucam_sanyal/metadata.csv"
sample_names_col <- 'Sample.name'
batches <- c("Sex", "Dataset")
template_df = read.csv(file = template_filename, sep=',')
rownames(template_df) = template_df[,sample_names_col]
template_df <- template_df[rownames(sorted_samples_df),]

# 3. Counts matrix. Sample names should be the same as in metadata file.
counts_filename <- "../../data/ucam_sanyal/counts_matrix.csv"
counts_matrix_df <- read.csv(file = counts_filename, sep = ",")
rownames(counts_matrix_df) = as.character(counts_matrix_df$X)
counts_matrix_df = counts_matrix_df[,-c(1)]
colnames(counts_matrix_df) <- gsub('\\.', ' ', colnames(counts_matrix_df))
counts_matrix_df <- counts_matrix_df[, rownames(sorted_samples_df)]
counts_matrix_df <- counts_matrix_df[(rowSums(counts_matrix_df)>dim(counts_matrix_df)[2]),] #Exclude low expressed counts
counts_matrix = data.matrix(counts_matrix_df)
colnames(counts_matrix) == row.names(template_df)

# 4. Output folder
output_dir <- '../../results/ucam_sanyal/finding_of_optimal_SWs_sequence/'


################################################################################
# Main process
################################################################################

# Parameters for the sw_graph
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
seq_of_sizes <- seq(18,20)
overlap <- 0.5


################################################################################
# 1. Build the sw_graph
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
save(sw_graph, file=paste(output_dir, "sw_graph_initial.RData", sep=''))



################################################################################
# 2. Scan the sw_graph by running the deseq worker for each edge
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
        degs <- deseq_worker(sw1_samples, sw2_samples, batches)
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
  filename = paste(output_dir, "sw_graph_tmp_", as.character(step-1), ".RData", sep='')
  save(sw_graph, runs, file=filename)
  
}

save(sw_graph, runs, file=paste(output_dir, "sw_graph_final.RData", sep=''))
