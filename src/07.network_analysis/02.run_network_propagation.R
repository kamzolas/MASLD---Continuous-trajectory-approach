source("../library.R")
library(igraph) # 2.2.1


################################################################################
# Description
################################################################################
# This is the main script to run network propagation for each sliding window
# and direction of differentiation (up & down). The user is able to define the 
# damping factor, the kind of normalization on network edges as well as how
# the null distribution is generated to perform statistical inference for the 
# propagation scores.
# Output files:
# - A new folder (propagation_results) which includes the results in RData format.
################################################################################


################################################################################
# A function to define the seed node weights for each sw.
################################################################################
define_personalized_values <- function(tmp_gene_scores, tmp_tf_scores, nodes) {
    # tmp_gene_scores: scores related to the results of the DE analysis
    # tmp_tf_scores: scores related to the results of the TF activity analysis
    # nodes: a list of network nodes to get and use the scores only for them
    genes_intersection <- setdiff(intersect(nodes, rownames(tmp_gene_scores)), rownames(tmp_tf_scores))
    tfs_intersection <- intersect(nodes, rownames(tmp_tf_scores))
    # Create the vector with length equal to network size
    personalized_values <- rep(1, length(nodes))
    names(personalized_values) <- nodes
    # Define the values for the input genes vector
    tmp_scores <- abs(tmp_gene_scores[genes_intersection,1])
    personalized_values[genes_intersection] <- tmp_scores
    # Define the values for the input tf vector
    tmp_scores <- abs(tmp_tf_scores[tfs_intersection,1])
    personalized_values[tfs_intersection] <- tmp_scores
    return(personalized_values)
}


################################################################################
# A function to define the random seed node scores. The score are derived from
# the 'distribution' parameter. A matrix of dimensionality N is returned, in 
# which each column corresponds to one random network propagation run.
################################################################################
define_random_personalized_values <- function(nodes, distribution, N=1000) {
    set.seed(1234)
    random_personalized_values_df <- sapply(1:N, function(n) {
        # get the random scores from the global distribution of scores for all
        # the SWs
        sample(x=distribution, size=length(nodes), replace=TRUE)
    })
    #random_noise_df <- sapply(1:N, function(n) {
      # get the random scores from the global distribution of scores for all
      # the SWs
    #  rtri(n=length(nodes), min=0, max=0.699, mode=0.523)
    #})
    #random_personalized_values_df <- random_personalized_values_df + random_noise_df
    rownames(random_personalized_values_df) <- nodes
    return(random_personalized_values_df)
}


################################################################################
# Personalized pagerank worker
################################################################################
pagerank_worker <- function(network, personalized, damping, directed) {
  pagerank_values <- page_rank(graph = network, 
                               directed = directed,
                               personalized = personalized, 
                               damping = damping,
                               weights = edge_attr(network, name=edge_weight_attr))
  pagerank_values <- pagerank_values$vector
  return(pagerank_values)
}


################################################################################
# Network propagation worker
################################################################################
network_propagation_worker <- function(scores, network, 
                                       values_distribution, damping=0.85, 
                                       directed=FALSE, N=1000) {

    ############################################################################
    # 1. Definition of network nodes
    ############################################################################
    nodes <- as_ids(V(network))
    
    ############################################################################
    # 2. Definition of the personalized values of network nodes
    ############################################################################
    # gene_scores and tf_scores contain the seed scores for all genes and TFs in 
    # the dataset. So the vector of network nodes need to be inserted as argument
    # in order to retrieve the seed scores only for this subset of dataset
    #personalized <- define_personalized_values(gene_scores, tf_scores, nodes)
    personalized <- scores[,1]
    names(personalized) <- rownames(scores)
    ############################################################################
    # 3. PageRank for the empirical personalized values
    ############################################################################
    empirical <- pagerank_worker(network = network,
                                 personalized = personalized, 
                                 damping = damping, 
                                 directed = directed)
    empirical_df <- data.frame(empirical)
    
    ############################################################################
    # 4.Construction of the matrix with randomized personalized values
    ############################################################################
    random_personalized_df <- define_random_personalized_values(nodes, values_distribution, N)
    
    ############################################################################
    # 5. PageRank for the randomized personalized values
    ############################################################################
    random_propagation_df = data.frame(x=nodes)
    rownames(random_propagation_df) <- nodes
    random_propagation_df$x <- NULL
    for (n in 1:N) {

      # Select column n with the randomized scores
      tmp_random_personalized <- random_personalized_df[,n]
      # Call of the PageRank worker
      tmp_random_propagation <- pagerank_worker(network = network,
                                                personalized = tmp_random_personalized, 
                                                damping = damping, 
                                                directed = directed)
      # Concatenation of the derived results with the main results matrix
      tmp_random_propagation <- as.data.frame(tmp_random_propagation)
      colnames(tmp_random_propagation) <- c(paste('random_', n, sep=''))
      random_propagation_df <- merge(random_propagation_df, tmp_random_propagation, 
                                     by=0, all=TRUE)
      rownames(random_propagation_df) <- random_propagation_df$Row.names
      random_propagation_df$Row.names <- NULL
    }
  
    ############################################################################
    # 6. Concatenation of empirical and random results
    ############################################################################
    output_df1 <- merge(random_propagation_df, empirical_df, by=0, all=TRUE)
    rownames(output_df1) <- output_df1$Row.names
    output_df1$Row.names <- NULL
    
    ############################################################################
    # 7. Calculation of p-values and return
    ############################################################################
    scores <- rowSums(output_df1[,dim(output_df1)[2]] > output_df1[,1:dim(output_df1)[2]-1])
    output_df2 <- data.frame(pvalue = (N-scores)/N)
    rownames(output_df2) <- rownames(output_df1)
    output_df2 <- output_df2[order(output_df2$pvalue, decreasing = FALSE),,drop=FALSE]
    output_df1 <- output_df1[,dim(output_df1)[2],drop=FALSE]
    return(list('propagation_values'=output_df1, 'pvalues'=output_df2))
}


################################################################################
# Inputs
################################################################################
#args = commandArgs(trailingOnly=TRUE)
#edge_weight_attr = args[1] #lapl_norm or lapl_norm_weight
#randomization_method = args[2] #sw or all
edge_weight_attr = 'lapl_norm_weight' #lapl_norm or lapl_norm_weight
randomization_method = 'sw' #sw or all
damping = 0.5
input_dir = '../../results/ucam_sanyal/networks/'
output_dir = '../../results/ucam_sanyal/networks/network_analysis/'


################################################################################
# 1. Loading of the network and initial weights (for DEA and TF activity 
# analysis).
################################################################################
load(paste(input_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- MASLD_unified_undirected_network

filename <- paste(input_dir, 'initial_weights/up.tsv', sep='')
up_scores <- read.csv(filename, sep='\t')

filename <- paste(input_dir, 'initial_weights/down.tsv', sep='')
down_scores <- read.csv(filename, sep='\t')

# The 'values_matrix' variable will be used for the randomization of weights.
values_matrix1 <- abs(as.matrix(up_scores))
values_matrix2 <- abs(as.matrix(down_scores))
values_matrix <- rbind(values_matrix1, values_matrix2)


################################################################################
# 2. Network analysis for each sw.
################################################################################
sw <- 'SW2'
SWs <- paste('SW', seq(2,dim(up_scores)[2]+1), sep='')
results <- list()
for (sw in SWs) {

    # Get the scores for this sw
    sw_up_scores <- up_scores[,sw,drop=FALSE]
    sw_down_scores <- down_scores[,sw,drop=FALSE]
    
    # Get the values for the randomization
    if (randomization_method == 'all') {
      values_distribution <- as.vector(values_matrix)
    } else if  (randomization_method == 'sw') {
      values_distribution <- as.vector(values_matrix[,sw])
    } else {
      print('error: no randomization method is defined')
    }

    # Run the propagation
    down_propagation <- network_propagation_worker(sw_down_scores,
                                                   network, values_distribution,
                                                   damping=damping, 
                                                   directed=FALSE, N=1000)

    up_propagation <- network_propagation_worker(sw_up_scores, 
                                                 network, values_distribution,
                                                 damping=damping,
                                                 directed=FALSE, N=1000)
    
    results[[sw]] <- list('up'=up_propagation, 'down'=down_propagation)
}


################################################################################
# 3. Save the results
################################################################################
propagation_results = results
output_dir = paste(output_dir, 'propagation_results/', sep='')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
filename <- paste(paste('random', randomization_method, sep='_'), 
                  edge_weight_attr, gsub("\\.", '', as.character(damping)), sep='-')
filename <- paste(output_dir, filename, '.RData', sep='')
save(propagation_results, file=filename)











