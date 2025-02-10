source("../library.R")
library(igraph) # ‘1.4.3’


################################################################################
#
# Output files:
#
# A new folder will be created ('network_analysis/initial_weights/') to store 
# two dataframes:
#
# - up.tsv : personalized 'upregulated' weights for network nodes in all sws 
# - down.tsv: personalized 'downregulated' weights for network nodes in all sws
#
################################################################################


define_gene_scores <- function(de_results) {
  SWs <- paste('SW', seq(2, (length(colnames(de_results))/2)+1), sep='')
  output_df <- data.frame(genes=rownames(de_results))
  rownames(output_df) <- rownames(de_results)
  sw <- 'SW2'
  for (sw in SWs) {
    col1 <- paste('logfc', sw, sep='_')
    col2 <- paste('fdr', sw, sep='_')
    tmp <- de_results[,c(col1, col2)]
    tmp[, col1] <- sign(tmp[, col1]) + tmp[, col1]
    tmp[, col2] <- -log10(tmp[, col2])
    tmp[, sw] <- tmp[, col1] * tmp[, col2]
    output_df <- merge(output_df, tmp[,sw,drop=FALSE], by.x=0, by.y=0, all=TRUE)
    rownames(output_df) <- output_df$Row.names
    output_df$Row.names <- NULL
  }
  output_df$genes <- NULL
  output_df[output_df > 10] <- 10
  output_df[output_df < -10] <- -10
  return(output_df)
}



define_tf_scores <- function(msviper_results) {
  SWs <- paste('SW', seq(2, (length(colnames(msviper_results))/2)+1), sep='')
  output_df <- data.frame(genes=rownames(msviper_results))
  rownames(output_df) <- rownames(msviper_results)
  sw <- 'SW2'
  for (sw in SWs) {
    col1 <- paste('nes', sw, sep='_')
    col2 <- paste('fdr', sw, sep='_')
    tmp <- msviper_results[,c(col1, col2)]
    tmp[, col1] <- sign(tmp[, col1]) 
    tmp[, col2] <- -log10(tmp[, col2])
    tmp[, sw] <- tmp[, col1] * tmp[, col2]
    output_df <- merge(output_df, tmp[,sw,drop=FALSE], by.x=0, by.y=0, all=TRUE)
    rownames(output_df) <- output_df$Row.names
    output_df$Row.names <- NULL
  }
  output_df$genes <- NULL
  #output_df[output_df > 10] <- 10
  #output_df[output_df < -10] <- -10
  return(output_df)
}





################################################################################
#
# A function to define the seed node weights for each sw.
#
################################################################################
define_personalized_values <- function(tmp_gene_scores, tmp_tf_scores, nodes) {
  
  # tmp_gene_scores: scores related to the results of the DE analysis
  # tmp_tf_scores: scores related to the results of the TF activity analysis
  # nodes: a list of network nodes to get and use the scores only for them
  genes_intersection <- setdiff(intersect(nodes, rownames(tmp_gene_scores)), 
                                rownames(tmp_tf_scores))
  tfs_intersection <- intersect(nodes, rownames(tmp_tf_scores))
  
  # Create a vector with length equal to network size
  personalized_values <- rep(1, length(nodes))
  #personalized_values <- rep(0, length(nodes))
  names(personalized_values) <- nodes
  
  # Define the values for the input genes
  tmp_scores <- 1 + abs(tmp_gene_scores[genes_intersection,1])
  personalized_values[genes_intersection] <- tmp_scores
  
  # Define the values for the input tfs
  tmp_scores <- 1 + abs(tmp_tf_scores[tfs_intersection,1])
  personalized_values[tfs_intersection] <- tmp_scores
  
  personalized_values[personalized_values > 10] <- 10
  
  return(personalized_values)
}



################################################################################
#
# 1. Input arguments
#
################################################################################
key_dir = '../../results/pseudo_bulk_analysis/networks/'



################################################################################
#
# 2. Loading of the unified MASLD network and the results from DEA and TF 
# activity analysis.
#
################################################################################
main_project_dir <- '../../results/networks/'
load(paste(main_project_dir, 'MASLD_unified_undirected_network_with_semantics.RData', sep=''))
network <- MASLD_unified_undirected_network

filename <- paste(key_dir, 'de_and_msviper/de_results.tsv', sep='')
de_results <- read.csv(filename, sep='\t')
gene_scores <- define_gene_scores(de_results)

filename <- paste(key_dir, 'de_and_msviper/msviper_results.tsv', sep='')
msviper_results <- read.csv(filename, sep='\t')
tf_scores <- define_tf_scores(msviper_results)


################################################################################
#
# 3. Creation of vectors with the personalized weights of network nodes, for 
# each sw and direction (up & down).
#
################################################################################
SWs <- paste('SW', seq(2,dim(tf_scores)[2]+1), sep='')
results <- list()
for (sw in SWs) {
  
  # Scores for the examined sw.
  sw_gene_scores <- gene_scores[,sw,drop=FALSE]
  sw_tf_scores <- tf_scores[,sw,drop=FALSE]
  
  # Separation the scores in positive and negative.
  down_genes_scores <- sw_gene_scores[sw_gene_scores[,sw] < 0,,drop=FALSE]
  down_tf_scores <- sw_tf_scores[,sw,drop=FALSE]
  #down_tf_scores[down_tf_scores > 0] <- -1
  down_tf_scores[down_tf_scores > 0] <- 0
  
  up_genes_scores <- sw_gene_scores[sw_gene_scores[,sw] > 0,,drop=FALSE]
  up_tf_scores <- sw_tf_scores[,sw,drop=FALSE]
  #up_tf_scores[up_tf_scores < 0] <- 1
  up_tf_scores[up_tf_scores < 0] <- 0
  
  # Calculation of the vectors with the personalized weights.
  nodes <- as_ids(V(network))
  down_personalized <- define_personalized_values(down_genes_scores, down_tf_scores, nodes)
  results[['down']][[sw]] <- down_personalized
  up_personalized <- define_personalized_values(up_genes_scores, up_tf_scores, nodes)
  results[['up']][[sw]] <- up_personalized
}

output_dir <- paste(key_dir, 'network_analysis/initial_weights/', sep='')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

filename <- paste(output_dir, 'up.tsv', sep='')
write.table(as.data.frame(results[['up']]), file=filename, quote=FALSE, sep='\t')

filename <- paste(output_dir, 'down.tsv', sep='')
write.table(as.data.frame(results[['down']]), file=filename, quote=FALSE, sep='\t')
